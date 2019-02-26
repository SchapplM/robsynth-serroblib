% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:12:38
% EndTime: 2019-02-26 22:12:38
% DurationCPUTime: 0.26s
% Computational Cost: add. (386->69), mult. (652->111), div. (0->0), fcn. (833->12), ass. (0->45)
t40 = sin(qJ(2));
t41 = sin(qJ(1));
t44 = cos(qJ(2));
t45 = cos(qJ(1));
t53 = cos(pkin(6));
t50 = t45 * t53;
t26 = t40 * t50 + t41 * t44;
t35 = qJ(3) + pkin(11);
t33 = sin(t35);
t34 = cos(t35);
t36 = sin(pkin(6));
t57 = t36 * t45;
t14 = t26 * t34 - t33 * t57;
t25 = t41 * t40 - t44 * t50;
t38 = sin(qJ(5));
t42 = cos(qJ(5));
t1 = t14 * t38 - t25 * t42;
t67 = t14 * t42 + t25 * t38;
t43 = cos(qJ(3));
t32 = t43 * pkin(3) + pkin(2);
t64 = pkin(10) + r_i_i_C(2);
t66 = pkin(4) * t34 + t64 * t33 + t32;
t54 = r_i_i_C(3) + qJ(6);
t65 = pkin(5) + r_i_i_C(1);
t46 = t54 * t38 + t65 * t42 + pkin(4);
t61 = t34 * t38;
t60 = t34 * t42;
t59 = t36 * t40;
t58 = t36 * t41;
t56 = t38 * t44;
t55 = t42 * t44;
t51 = t41 * t53;
t39 = sin(qJ(3));
t49 = t36 * (pkin(3) * t39 + pkin(8));
t48 = -t26 * t33 - t34 * t57;
t37 = -qJ(4) - pkin(9);
t28 = -t40 * t51 + t45 * t44;
t27 = t45 * t40 + t44 * t51;
t22 = t53 * t33 + t34 * t59;
t18 = t28 * t34 + t33 * t58;
t17 = t28 * t33 - t34 * t58;
t11 = t22 * t38 + t36 * t55;
t6 = t18 * t42 + t27 * t38;
t5 = t18 * t38 - t27 * t42;
t2 = [-t41 * pkin(1) - t14 * pkin(4) - t54 * t1 + t25 * t37 - t26 * t32 + t45 * t49 + t64 * t48 - t65 * t67, -t28 * t37 + t54 * (-t27 * t61 - t28 * t42) + t65 * (-t27 * t60 + t28 * t38) - t66 * t27, t64 * t18 + (-t28 * t39 + t43 * t58) * pkin(3) - t46 * t17, t27, -t65 * t5 + t54 * t6, t5; t45 * pkin(1) + t18 * pkin(4) + t64 * t17 - t27 * t37 + t28 * t32 + t41 * t49 + t54 * t5 + t65 * t6, -t26 * t37 + t65 * (-t25 * t60 + t26 * t38) + t54 * (-t25 * t61 - t26 * t42) - t66 * t25, t64 * t14 + (-t26 * t39 - t43 * t57) * pkin(3) + t46 * t48, t25, -t65 * t1 + t54 * t67, t1; 0 (t65 * (t34 * t55 + t38 * t40) + t54 * (t34 * t56 - t40 * t42) - t37 * t40 + t66 * t44) * t36, t64 * t22 + (-t39 * t59 + t53 * t43) * pkin(3) + t46 * (-t33 * t59 + t53 * t34) -t36 * t44, t54 * (t22 * t42 - t36 * t56) - t65 * t11, t11;];
Ja_transl  = t2;
