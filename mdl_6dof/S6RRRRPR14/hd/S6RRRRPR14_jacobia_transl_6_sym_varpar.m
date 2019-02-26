% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR14
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR14_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR14_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:38:11
% EndTime: 2019-02-26 22:38:11
% DurationCPUTime: 0.35s
% Computational Cost: add. (658->92), mult. (1640->161), div. (0->0), fcn. (2142->16), ass. (0->62)
t55 = sin(qJ(2));
t57 = cos(qJ(2));
t58 = cos(qJ(1));
t79 = cos(pkin(6));
t71 = t58 * t79;
t85 = sin(qJ(1));
t37 = t85 * t55 - t57 * t71;
t38 = t55 * t71 + t85 * t57;
t54 = sin(qJ(3));
t78 = cos(pkin(7));
t72 = t54 * t78;
t50 = sin(pkin(7));
t51 = sin(pkin(6));
t81 = t51 * t58;
t76 = t50 * t81;
t86 = cos(qJ(3));
t16 = -t37 * t72 + t38 * t86 - t54 * t76;
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t73 = t51 * t78;
t65 = t37 * t50 - t58 * t73;
t4 = t16 * t56 + t65 * t53;
t3 = t16 * t53 - t65 * t56;
t68 = t79 * t85;
t63 = t58 * t55 + t57 * t68;
t74 = t51 * t85;
t89 = -t50 * t74 + t63 * t78;
t88 = pkin(10) * t50;
t87 = r_i_i_C(3) + pkin(12) + qJ(5);
t84 = t50 * t53;
t83 = t50 * t56;
t82 = t51 * t57;
t80 = t55 * t50;
t77 = t51 * t80;
t75 = sin(pkin(13)) * pkin(5) + pkin(11);
t70 = t79 * t50;
t67 = t78 * t86;
t45 = cos(pkin(13)) * pkin(5) + pkin(4);
t48 = pkin(13) + qJ(6);
t46 = sin(t48);
t47 = cos(t48);
t66 = t47 * r_i_i_C(1) - t46 * r_i_i_C(2) + t45;
t64 = t46 * r_i_i_C(1) + t47 * r_i_i_C(2) + t75;
t62 = -t50 * t82 + t79 * t78;
t60 = -t87 * t53 - t66 * t56 - pkin(3);
t15 = t37 * t67 + t38 * t54 + t86 * t76;
t59 = t63 * t50 + t85 * t73;
t39 = -t55 * t68 + t58 * t57;
t35 = (-t55 * t72 + t86 * t57) * t51;
t30 = t54 * t70 + (t86 * t55 + t57 * t72) * t51;
t29 = t51 * t55 * t54 - t67 * t82 - t86 * t70;
t24 = -t39 * t72 - t63 * t86;
t22 = -t37 * t86 - t38 * t72;
t20 = t39 * t86 - t89 * t54;
t19 = t39 * t54 + t89 * t86;
t14 = t30 * t56 + t62 * t53;
t13 = t30 * t53 - t62 * t56;
t8 = t20 * t56 + t59 * t53;
t7 = t20 * t53 - t59 * t56;
t2 = t19 * t46 + t8 * t47;
t1 = t19 * t47 - t8 * t46;
t5 = [-t85 * pkin(1) - t38 * pkin(2) - t16 * pkin(3) + pkin(9) * t81 - t65 * pkin(10) - t64 * t15 - t87 * t3 - t4 * t66, t24 * pkin(3) - t63 * pkin(2) + t39 * t88 + t66 * (t24 * t56 + t39 * t84) + t64 * (t39 * t67 - t63 * t54) + t87 * (t24 * t53 - t39 * t83) t60 * t19 + t64 * t20, -t66 * t7 + t87 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t58 * pkin(1) + t39 * pkin(2) + t20 * pkin(3) + pkin(9) * t74 + t59 * pkin(10) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t75 * t19 + t8 * t45 + t87 * t7, t38 * t88 - t37 * pkin(2) + t22 * pkin(3) + t87 * (t22 * t53 - t38 * t83) + t66 * (t22 * t56 + t38 * t84) + t64 * (-t37 * t54 + t38 * t67) t60 * t15 + t64 * t16, -t66 * t3 + t87 * t4, t3 (t15 * t47 - t4 * t46) * r_i_i_C(1) + (-t15 * t46 - t4 * t47) * r_i_i_C(2); 0, t35 * pkin(3) + t66 * (t35 * t56 + t53 * t77) + t87 * (t35 * t53 - t56 * t77) + (t57 * pkin(2) + pkin(10) * t80 + t64 * (t54 * t57 + t55 * t67)) * t51, t60 * t29 + t64 * t30, -t66 * t13 + t87 * t14, t13 (-t14 * t46 + t29 * t47) * r_i_i_C(1) + (-t14 * t47 - t29 * t46) * r_i_i_C(2);];
Ja_transl  = t5;
