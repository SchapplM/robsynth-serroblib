% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP1
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
% Datum: 2019-02-26 19:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPRRP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:50:18
% EndTime: 2019-02-26 19:50:18
% DurationCPUTime: 0.21s
% Computational Cost: add. (253->41), mult. (637->75), div. (0->0), fcn. (836->12), ass. (0->37)
t54 = pkin(5) + r_i_i_C(1);
t25 = sin(pkin(11));
t33 = sin(qJ(2));
t36 = cos(qJ(2));
t48 = cos(pkin(11));
t40 = -t33 * t25 + t36 * t48;
t31 = sin(qJ(5));
t34 = cos(qJ(5));
t55 = t34 * r_i_i_C(2) + t54 * t31 + pkin(8);
t32 = sin(qJ(4));
t35 = cos(qJ(4));
t41 = -r_i_i_C(2) * t31 + t54 * t34 + pkin(4);
t53 = r_i_i_C(3) + qJ(6) + pkin(9);
t37 = t53 * t32 + t41 * t35 + pkin(3);
t26 = sin(pkin(10));
t27 = sin(pkin(6));
t52 = t26 * t27;
t28 = cos(pkin(10));
t51 = t28 * t27;
t29 = cos(pkin(6));
t50 = t29 * t36;
t19 = -t36 * t25 - t33 * t48;
t17 = t19 * t29;
t7 = -t17 * t28 + t26 * t40;
t42 = -t26 * t17 - t28 * t40;
t38 = t29 * t40;
t16 = t19 * t27;
t15 = t40 * t27;
t12 = -t16 * t35 + t29 * t32;
t11 = -t16 * t32 - t29 * t35;
t9 = t19 * t28 - t26 * t38;
t6 = t26 * t19 + t28 * t38;
t4 = t32 * t52 - t35 * t42;
t3 = -t32 * t42 - t35 * t52;
t2 = -t32 * t51 + t35 * t7;
t1 = t32 * t7 + t35 * t51;
t5 = [0 (-t26 * t50 - t28 * t33) * pkin(2) - t55 * t42 + t37 * t9, t52, -t41 * t3 + t53 * t4 (t31 * t9 - t34 * t4) * r_i_i_C(2) + t54 * (-t31 * t4 - t34 * t9) t3; 0 (-t26 * t33 + t28 * t50) * pkin(2) + t55 * t7 + t37 * t6, -t51, -t41 * t1 + t53 * t2 (-t2 * t34 + t31 * t6) * r_i_i_C(2) + t54 * (-t2 * t31 - t34 * t6) t1; 1, t27 * t36 * pkin(2) + t37 * t15 - t55 * t16, t29, -t41 * t11 + t53 * t12 (-t12 * t34 + t15 * t31) * r_i_i_C(2) + t54 * (-t12 * t31 - t15 * t34) t11;];
Ja_transl  = t5;
