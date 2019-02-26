% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:04:05
% EndTime: 2019-02-26 20:04:06
% DurationCPUTime: 0.19s
% Computational Cost: add. (351->46), mult. (405->71), div. (0->0), fcn. (506->14), ass. (0->37)
t59 = pkin(10) + r_i_i_C(3);
t35 = sin(qJ(6));
t37 = cos(qJ(6));
t58 = t37 * r_i_i_C(1) - t35 * r_i_i_C(2) + pkin(5);
t31 = qJ(3) + pkin(12);
t23 = pkin(4) * cos(t31) + cos(qJ(3)) * pkin(3);
t28 = qJ(5) + t31;
t26 = sin(t28);
t27 = cos(t28);
t57 = t59 * t26 + t58 * t27 + pkin(2) + t23;
t32 = sin(pkin(11));
t33 = sin(pkin(6));
t53 = t32 * t33;
t36 = sin(qJ(2));
t52 = t32 * t36;
t38 = cos(qJ(2));
t51 = t32 * t38;
t50 = t33 * t36;
t49 = t33 * t38;
t48 = cos(pkin(11));
t47 = t33 * t48;
t46 = t48 * t36;
t45 = t48 * t38;
t43 = t35 * r_i_i_C(1) + t37 * r_i_i_C(2) + pkin(8) + pkin(9) + qJ(4);
t34 = cos(pkin(6));
t17 = t34 * t46 + t51;
t8 = t17 * t27 - t26 * t47;
t42 = t59 * t8 + t58 * (-t17 * t26 - t27 * t47);
t19 = -t34 * t52 + t45;
t10 = t19 * t27 + t26 * t53;
t41 = t59 * t10 + t58 * (-t19 * t26 + t27 * t53);
t15 = t34 * t26 + t27 * t50;
t40 = t59 * t15 + t58 * (-t26 * t50 + t34 * t27);
t22 = -pkin(4) * sin(t31) - sin(qJ(3)) * pkin(3);
t18 = t34 * t51 + t46;
t16 = -t34 * t45 + t52;
t1 = [0, -t18 * t57 + t43 * t19, t19 * t22 + t23 * t53 + t41, t18, t41 (-t10 * t35 + t18 * t37) * r_i_i_C(1) + (-t10 * t37 - t18 * t35) * r_i_i_C(2); 0, -t16 * t57 + t43 * t17, t17 * t22 - t23 * t47 + t42, t16, t42 (t16 * t37 - t8 * t35) * r_i_i_C(1) + (-t16 * t35 - t8 * t37) * r_i_i_C(2); 1 (t43 * t36 + t57 * t38) * t33, t22 * t50 + t34 * t23 + t40, -t49, t40 (-t15 * t35 - t37 * t49) * r_i_i_C(1) + (-t15 * t37 + t35 * t49) * r_i_i_C(2);];
Ja_transl  = t1;
