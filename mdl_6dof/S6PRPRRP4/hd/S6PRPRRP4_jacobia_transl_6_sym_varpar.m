% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRP4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:52:00
% EndTime: 2019-02-26 19:52:01
% DurationCPUTime: 0.19s
% Computational Cost: add. (276->48), mult. (463->84), div. (0->0), fcn. (594->11), ass. (0->37)
t29 = pkin(11) + qJ(4);
t27 = sin(t29);
t28 = cos(t29);
t52 = pkin(9) + r_i_i_C(2);
t54 = t52 * t27 + pkin(4) * t28 + cos(pkin(11)) * pkin(3) + pkin(2);
t53 = pkin(5) + r_i_i_C(1);
t33 = sin(qJ(5));
t51 = t28 * t33;
t35 = cos(qJ(5));
t50 = t28 * t35;
t30 = sin(pkin(10));
t31 = sin(pkin(6));
t49 = t30 * t31;
t34 = sin(qJ(2));
t48 = t31 * t34;
t36 = cos(qJ(2));
t47 = t33 * t36;
t46 = t35 * t36;
t45 = r_i_i_C(3) + qJ(6);
t44 = cos(pkin(6));
t43 = cos(pkin(10));
t41 = t30 * t44;
t40 = t31 * t43;
t39 = t44 * t43;
t37 = t45 * t33 + t53 * t35 + pkin(4);
t32 = -pkin(8) - qJ(3);
t24 = -t34 * t41 + t43 * t36;
t23 = t43 * t34 + t36 * t41;
t22 = t30 * t36 + t34 * t39;
t21 = t30 * t34 - t36 * t39;
t18 = t44 * t27 + t28 * t48;
t13 = t18 * t33 + t31 * t46;
t12 = t24 * t28 + t27 * t49;
t10 = t22 * t28 - t27 * t40;
t3 = t12 * t33 - t23 * t35;
t1 = t10 * t33 - t21 * t35;
t2 = [0, -t24 * t32 + t53 * (-t23 * t50 + t24 * t33) + t45 * (-t23 * t51 - t24 * t35) - t54 * t23, t23, t52 * t12 + t37 * (-t24 * t27 + t28 * t49) t45 * (t12 * t35 + t23 * t33) - t53 * t3, t3; 0, -t22 * t32 + t53 * (-t21 * t50 + t22 * t33) + t45 * (-t21 * t51 - t22 * t35) - t54 * t21, t21, t52 * t10 + t37 * (-t22 * t27 - t28 * t40) t45 * (t10 * t35 + t21 * t33) - t53 * t1, t1; 1 (t53 * (t28 * t46 + t33 * t34) + t45 * (t28 * t47 - t34 * t35) - t32 * t34 + t54 * t36) * t31, -t31 * t36, t52 * t18 + t37 * (-t27 * t48 + t44 * t28) t45 * (t18 * t35 - t31 * t47) - t53 * t13, t13;];
Ja_transl  = t2;
