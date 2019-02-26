% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_jacobia_transl_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:10:38
% EndTime: 2019-02-26 20:10:38
% DurationCPUTime: 0.20s
% Computational Cost: add. (365->52), mult. (426->82), div. (0->0), fcn. (529->14), ass. (0->37)
t59 = pkin(10) + r_i_i_C(3);
t38 = sin(qJ(6));
t40 = cos(qJ(6));
t58 = t40 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(5);
t33 = qJ(3) + qJ(4);
t30 = cos(t33);
t23 = pkin(4) * t30 + cos(qJ(3)) * pkin(3);
t28 = pkin(12) + t33;
t25 = sin(t28);
t26 = cos(t28);
t57 = t59 * t25 + t58 * t26 + pkin(2) + t23;
t34 = sin(pkin(11));
t35 = sin(pkin(6));
t53 = t34 * t35;
t36 = cos(pkin(11));
t52 = t35 * t36;
t39 = sin(qJ(2));
t51 = t35 * t39;
t41 = cos(qJ(2));
t50 = t35 * t41;
t37 = cos(pkin(6));
t49 = t37 * t39;
t48 = t37 * t41;
t46 = t38 * r_i_i_C(1) + t40 * r_i_i_C(2) + pkin(8) + pkin(9) + qJ(5);
t17 = t34 * t41 + t36 * t49;
t8 = t17 * t26 - t25 * t52;
t45 = t59 * t8 + t58 * (-t17 * t25 - t26 * t52);
t19 = -t34 * t49 + t36 * t41;
t10 = t19 * t26 + t25 * t53;
t44 = t59 * t10 + t58 * (-t19 * t25 + t26 * t53);
t15 = t37 * t25 + t26 * t51;
t43 = t59 * t15 + t58 * (-t25 * t51 + t37 * t26);
t29 = sin(t33);
t22 = -pkin(4) * t29 - sin(qJ(3)) * pkin(3);
t18 = t34 * t48 + t36 * t39;
t16 = t34 * t39 - t36 * t48;
t1 = [0, -t18 * t57 + t46 * t19, t19 * t22 + t23 * t53 + t44 (-t19 * t29 + t30 * t53) * pkin(4) + t44, t18 (-t10 * t38 + t18 * t40) * r_i_i_C(1) + (-t10 * t40 - t18 * t38) * r_i_i_C(2); 0, -t16 * t57 + t46 * t17, t17 * t22 - t23 * t52 + t45 (-t17 * t29 - t30 * t52) * pkin(4) + t45, t16 (t16 * t40 - t8 * t38) * r_i_i_C(1) + (-t16 * t38 - t8 * t40) * r_i_i_C(2); 1 (t46 * t39 + t57 * t41) * t35, t22 * t51 + t37 * t23 + t43 (-t29 * t51 + t30 * t37) * pkin(4) + t43, -t50 (-t15 * t38 - t40 * t50) * r_i_i_C(1) + (-t15 * t40 + t38 * t50) * r_i_i_C(2);];
Ja_transl  = t1;
