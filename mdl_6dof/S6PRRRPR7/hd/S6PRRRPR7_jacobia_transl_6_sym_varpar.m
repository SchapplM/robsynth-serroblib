% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRRPR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:14:03
% EndTime: 2019-02-26 20:14:04
% DurationCPUTime: 0.29s
% Computational Cost: add. (504->75), mult. (1262->140), div. (0->0), fcn. (1646->16), ass. (0->57)
t41 = sin(pkin(7));
t47 = sin(qJ(2));
t49 = cos(qJ(2));
t43 = cos(pkin(12));
t71 = cos(pkin(6));
t65 = t43 * t71;
t69 = sin(pkin(12));
t56 = t69 * t47 - t49 * t65;
t70 = cos(pkin(7));
t42 = sin(pkin(6));
t74 = t42 * t43;
t81 = t41 * t74 + t56 * t70;
t60 = t71 * t69;
t57 = t43 * t47 + t49 * t60;
t66 = t42 * t69;
t80 = -t41 * t66 + t57 * t70;
t79 = pkin(9) * t41;
t78 = r_i_i_C(3) + pkin(11) + qJ(5);
t77 = cos(qJ(3));
t45 = sin(qJ(4));
t76 = t41 * t45;
t48 = cos(qJ(4));
t75 = t41 * t48;
t73 = t42 * t47;
t72 = t42 * t49;
t68 = t41 * t73;
t46 = sin(qJ(3));
t64 = t46 * t70;
t63 = t71 * t41;
t61 = t70 * t77;
t39 = pkin(13) + qJ(6);
t37 = sin(t39);
t38 = cos(t39);
t59 = t38 * r_i_i_C(1) - t37 * r_i_i_C(2) + cos(pkin(13)) * pkin(5) + pkin(4);
t58 = sin(pkin(13)) * pkin(5) + t37 * r_i_i_C(1) + t38 * r_i_i_C(2) + pkin(10);
t55 = -t41 * t72 + t71 * t70;
t54 = -t78 * t45 - t59 * t48 - pkin(3);
t51 = t56 * t41 - t70 * t74;
t50 = t57 * t41 + t70 * t66;
t31 = t43 * t49 - t47 * t60;
t30 = t47 * t65 + t69 * t49;
t28 = (-t47 * t64 + t77 * t49) * t42;
t24 = t46 * t63 + (t77 * t47 + t49 * t64) * t42;
t23 = t46 * t73 - t61 * t72 - t77 * t63;
t18 = -t31 * t64 - t57 * t77;
t16 = -t30 * t64 - t56 * t77;
t14 = t24 * t48 + t55 * t45;
t13 = t24 * t45 - t55 * t48;
t12 = t31 * t77 - t80 * t46;
t11 = t31 * t46 + t80 * t77;
t10 = t30 * t77 - t81 * t46;
t9 = t30 * t46 + t81 * t77;
t4 = t12 * t48 + t50 * t45;
t3 = t12 * t45 - t50 * t48;
t2 = t10 * t48 + t51 * t45;
t1 = t10 * t45 - t51 * t48;
t5 = [0, t18 * pkin(3) - t57 * pkin(2) + t31 * t79 + t78 * (t18 * t45 - t31 * t75) + t59 * (t18 * t48 + t31 * t76) + t58 * (t31 * t61 - t57 * t46) t54 * t11 + t58 * t12, -t59 * t3 + t78 * t4, t3 (t11 * t38 - t4 * t37) * r_i_i_C(1) + (-t11 * t37 - t4 * t38) * r_i_i_C(2); 0, t16 * pkin(3) - t56 * pkin(2) + t30 * t79 + t78 * (t16 * t45 - t30 * t75) + t59 * (t16 * t48 + t30 * t76) + t58 * (t30 * t61 - t56 * t46) t58 * t10 + t54 * t9, -t59 * t1 + t78 * t2, t1 (-t2 * t37 + t9 * t38) * r_i_i_C(1) + (-t2 * t38 - t9 * t37) * r_i_i_C(2); 1, t28 * pkin(3) + t59 * (t28 * t48 + t45 * t68) + t78 * (t28 * t45 - t48 * t68) + (t49 * pkin(2) + t47 * t79 + t58 * (t46 * t49 + t47 * t61)) * t42, t54 * t23 + t58 * t24, -t59 * t13 + t78 * t14, t13 (-t14 * t37 + t23 * t38) * r_i_i_C(1) + (-t14 * t38 - t23 * t37) * r_i_i_C(2);];
Ja_transl  = t5;
