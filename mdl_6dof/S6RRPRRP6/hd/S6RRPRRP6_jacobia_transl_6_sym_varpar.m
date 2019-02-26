% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP6_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP6_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:49
% EndTime: 2019-02-26 21:48:50
% DurationCPUTime: 0.27s
% Computational Cost: add. (456->69), mult. (1155->112), div. (0->0), fcn. (1534->12), ass. (0->46)
t45 = sin(pkin(11));
t50 = sin(qJ(2));
t54 = cos(qJ(2));
t62 = cos(pkin(11));
t39 = -t54 * t45 - t50 * t62;
t47 = cos(pkin(6));
t36 = t39 * t47;
t51 = sin(qJ(1));
t55 = cos(qJ(1));
t58 = -t50 * t45 + t54 * t62;
t25 = -t55 * t36 + t51 * t58;
t49 = sin(qJ(4));
t53 = cos(qJ(4));
t46 = sin(pkin(6));
t64 = t55 * t46;
t14 = t25 * t53 - t49 * t64;
t35 = t58 * t47;
t24 = t55 * t35 + t51 * t39;
t48 = sin(qJ(5));
t52 = cos(qJ(5));
t1 = t14 * t48 + t24 * t52;
t75 = t14 * t52 - t24 * t48;
t73 = pkin(10) + r_i_i_C(2);
t57 = pkin(4) * t53 + t73 * t49 + pkin(3);
t63 = r_i_i_C(3) + qJ(6);
t74 = pkin(5) + r_i_i_C(1);
t56 = t63 * t48 + t74 * t52 + pkin(4);
t72 = t54 * pkin(2);
t69 = t47 * t54;
t68 = t48 * t53;
t66 = t51 * t46;
t65 = t52 * t53;
t26 = -t51 * t36 - t55 * t58;
t59 = -t25 * t49 - t53 * t64;
t44 = pkin(1) + t72;
t37 = t47 * t50 * pkin(2) + (-pkin(8) - qJ(3)) * t46;
t34 = t39 * t46;
t33 = t58 * t46;
t30 = -t34 * t53 + t47 * t49;
t27 = -t51 * t35 + t55 * t39;
t18 = -t26 * t53 + t49 * t66;
t17 = -t26 * t49 - t53 * t66;
t11 = t30 * t48 + t33 * t52;
t6 = t18 * t52 - t27 * t48;
t5 = t18 * t48 + t27 * t52;
t2 = [-t25 * pkin(3) - t14 * pkin(4) + t24 * pkin(9) - t63 * t1 - t55 * t37 - t51 * t44 + t73 * t59 - t74 * t75, -t26 * pkin(9) + t63 * (t26 * t52 + t27 * t68) + t74 * (-t26 * t48 + t27 * t65) + (-t50 * t55 - t51 * t69) * pkin(2) + t57 * t27, t66, -t56 * t17 + t73 * t18, -t74 * t5 + t63 * t6, t5; -pkin(3) * t26 + t18 * pkin(4) - t27 * pkin(9) + t73 * t17 - t51 * t37 + t55 * t44 + t63 * t5 + t74 * t6, t25 * pkin(9) + t74 * (t24 * t65 + t25 * t48) + t63 * (t24 * t68 - t25 * t52) + (-t50 * t51 + t55 * t69) * pkin(2) + t57 * t24, -t64, t73 * t14 + t56 * t59, -t74 * t1 + t63 * t75, t1; 0, t46 * t72 - t34 * pkin(9) + t74 * (t33 * t65 - t34 * t48) + t63 * (t33 * t68 + t34 * t52) + t57 * t33, t47, t73 * t30 + t56 * (t34 * t49 + t47 * t53) t63 * (t30 * t52 - t33 * t48) - t74 * t11, t11;];
Ja_transl  = t2;
