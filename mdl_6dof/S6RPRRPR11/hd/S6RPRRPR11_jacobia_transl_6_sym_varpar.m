% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR11_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR11_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobia_transl_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:44
% EndTime: 2019-02-26 21:06:44
% DurationCPUTime: 0.31s
% Computational Cost: add. (509->66), mult. (1267->110), div. (0->0), fcn. (1667->16), ass. (0->53)
t42 = cos(qJ(1));
t69 = cos(pkin(12));
t71 = cos(pkin(6));
t59 = t71 * t69;
t67 = sin(pkin(12));
t73 = sin(qJ(1));
t53 = -t42 * t59 + t73 * t67;
t37 = sin(pkin(6));
t68 = sin(pkin(7));
t63 = t37 * t68;
t70 = cos(pkin(7));
t82 = t42 * t63 + t53 * t70;
t57 = t71 * t67;
t25 = t42 * t57 + t73 * t69;
t40 = sin(qJ(3));
t74 = cos(qJ(3));
t14 = -t25 * t74 + t82 * t40;
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t64 = t37 * t70;
t44 = -t42 * t64 + t53 * t68;
t81 = t14 * t39 + t44 * t41;
t80 = t14 * t41 - t44 * t39;
t11 = t25 * t40 + t82 * t74;
t48 = t42 * t67 + t73 * t59;
t76 = t48 * t70 - t73 * t63;
t75 = r_i_i_C(3) + pkin(11) + qJ(5);
t72 = t42 * t37;
t66 = sin(pkin(13)) * pkin(5) + pkin(10);
t65 = t73 * t37;
t58 = t71 * t68;
t56 = t70 * t69;
t32 = cos(pkin(13)) * pkin(5) + pkin(4);
t35 = pkin(13) + qJ(6);
t33 = sin(t35);
t34 = cos(t35);
t55 = t34 * r_i_i_C(1) - t33 * r_i_i_C(2) + t32;
t52 = t33 * r_i_i_C(1) + t34 * r_i_i_C(2) + t66;
t47 = -t75 * t39 - t55 * t41 - pkin(3);
t46 = -t69 * t63 + t71 * t70;
t43 = t48 * t68 + t73 * t64;
t26 = t42 * t69 - t73 * t57;
t20 = t40 * t58 + (t40 * t56 + t67 * t74) * t37;
t19 = -t74 * t58 + (t40 * t67 - t56 * t74) * t37;
t16 = t26 * t74 - t76 * t40;
t15 = t26 * t40 + t76 * t74;
t10 = t20 * t41 + t46 * t39;
t9 = t20 * t39 - t46 * t41;
t8 = t16 * t41 + t43 * t39;
t7 = t16 * t39 - t43 * t41;
t2 = t15 * t33 + t8 * t34;
t1 = t15 * t34 - t8 * t33;
t3 = [-t73 * pkin(1) - t25 * pkin(2) + t14 * pkin(3) - t44 * pkin(9) + qJ(2) * t72 - t52 * t11 + t55 * t80 + t75 * t81, t65, t47 * t15 + t52 * t16, -t55 * t7 + t75 * t8, t7, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t42 * pkin(1) + t26 * pkin(2) + t16 * pkin(3) + t43 * pkin(9) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + qJ(2) * t65 + t66 * t15 + t8 * t32 + t75 * t7, -t72, t47 * t11 - t14 * t52, t55 * t81 - t75 * t80, -t81 (t11 * t34 + t33 * t80) * r_i_i_C(1) + (-t11 * t33 + t34 * t80) * r_i_i_C(2); 0, t71, t47 * t19 + t52 * t20, t75 * t10 - t55 * t9, t9 (-t10 * t33 + t19 * t34) * r_i_i_C(1) + (-t10 * t34 - t19 * t33) * r_i_i_C(2);];
Ja_transl  = t3;
