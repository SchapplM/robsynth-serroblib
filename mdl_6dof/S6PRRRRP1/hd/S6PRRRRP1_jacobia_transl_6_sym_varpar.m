% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:15
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRRP1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRP1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:15:12
% EndTime: 2019-02-26 20:15:13
% DurationCPUTime: 0.23s
% Computational Cost: add. (293->42), mult. (460->74), div. (0->0), fcn. (573->12), ass. (0->39)
t63 = pkin(5) + r_i_i_C(1);
t65 = r_i_i_C(3) + qJ(6) + pkin(10);
t35 = sin(qJ(5));
t38 = cos(qJ(5));
t64 = -t35 * r_i_i_C(2) + t63 * t38 + pkin(4);
t30 = qJ(3) + qJ(4);
t28 = sin(t30);
t29 = cos(t30);
t39 = cos(qJ(3));
t62 = t39 * pkin(3) + t65 * t28 + t64 * t29 + pkin(2);
t31 = sin(pkin(11));
t32 = sin(pkin(6));
t58 = t31 * t32;
t33 = cos(pkin(11));
t57 = t32 * t33;
t37 = sin(qJ(2));
t56 = t32 * t37;
t55 = t32 * t39;
t40 = cos(qJ(2));
t54 = t32 * t40;
t53 = cos(pkin(6));
t52 = t37 * t53;
t51 = t40 * t53;
t20 = t31 * t40 + t33 * t52;
t10 = t20 * t29 - t28 * t57;
t9 = t20 * t28 + t29 * t57;
t46 = t65 * t10 - t64 * t9;
t22 = -t31 * t52 + t33 * t40;
t11 = t22 * t28 - t29 * t58;
t12 = t22 * t29 + t28 * t58;
t45 = -t64 * t11 + t65 * t12;
t44 = t38 * r_i_i_C(2) + t63 * t35 + pkin(8) + pkin(9);
t17 = t28 * t56 - t53 * t29;
t18 = t53 * t28 + t29 * t56;
t43 = -t64 * t17 + t65 * t18;
t36 = sin(qJ(3));
t21 = t31 * t51 + t33 * t37;
t19 = t31 * t37 - t33 * t51;
t1 = [0, -t21 * t62 + t44 * t22 (-t22 * t36 + t31 * t55) * pkin(3) + t45, t45 (-t12 * t38 - t21 * t35) * r_i_i_C(2) + t63 * (-t12 * t35 + t21 * t38) t11; 0, -t19 * t62 + t44 * t20 (-t20 * t36 - t33 * t55) * pkin(3) + t46, t46 (-t10 * t38 - t19 * t35) * r_i_i_C(2) + t63 * (-t10 * t35 + t19 * t38) t9; 1 (t44 * t37 + t62 * t40) * t32 (-t36 * t56 + t53 * t39) * pkin(3) + t43, t43 (-t18 * t38 + t35 * t54) * r_i_i_C(2) + t63 * (-t18 * t35 - t38 * t54) t17;];
Ja_transl  = t1;
