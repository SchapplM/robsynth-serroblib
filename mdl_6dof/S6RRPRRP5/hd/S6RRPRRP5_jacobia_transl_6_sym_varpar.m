% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP5
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
% Datum: 2019-02-26 21:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRRP5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:20
% EndTime: 2019-02-26 21:48:21
% DurationCPUTime: 0.25s
% Computational Cost: add. (355->58), mult. (869->94), div. (0->0), fcn. (1142->12), ass. (0->43)
t62 = pkin(5) + r_i_i_C(1);
t32 = sin(pkin(11));
t38 = sin(qJ(2));
t42 = cos(qJ(2));
t55 = cos(pkin(11));
t48 = -t38 * t32 + t42 * t55;
t36 = sin(qJ(5));
t40 = cos(qJ(5));
t46 = t40 * r_i_i_C(2) + t62 * t36 + pkin(9);
t37 = sin(qJ(4));
t41 = cos(qJ(4));
t30 = t40 * pkin(5) + pkin(4);
t49 = t40 * r_i_i_C(1) - t36 * r_i_i_C(2) + t30;
t60 = r_i_i_C(3) + qJ(6) + pkin(10);
t44 = t60 * t37 + t49 * t41 + pkin(3);
t61 = t42 * pkin(2);
t34 = cos(pkin(6));
t59 = t34 * t42;
t33 = sin(pkin(6));
t39 = sin(qJ(1));
t57 = t39 * t33;
t43 = cos(qJ(1));
t56 = t43 * t33;
t24 = -t42 * t32 - t38 * t55;
t21 = t24 * t34;
t11 = -t43 * t21 + t39 * t48;
t4 = t11 * t41 - t37 * t56;
t45 = t48 * t34;
t13 = t43 * t24 - t39 * t45;
t50 = -t39 * t21 - t43 * t48;
t8 = t37 * t57 - t41 * t50;
t1 = -t13 * t40 - t8 * t36;
t3 = t11 * t37 + t41 * t56;
t31 = pkin(1) + t61;
t22 = t34 * t38 * pkin(2) + (-pkin(8) - qJ(3)) * t33;
t20 = t24 * t33;
t19 = t48 * t33;
t16 = -t20 * t41 + t34 * t37;
t15 = -t20 * t37 - t34 * t41;
t10 = t39 * t24 + t43 * t45;
t7 = -t37 * t50 - t41 * t57;
t2 = -t13 * t36 + t8 * t40;
t5 = [-t11 * pkin(3) + t46 * t10 - t43 * t22 - t60 * t3 - t39 * t31 - t49 * t4 (-t43 * t38 - t39 * t59) * pkin(2) - t46 * t50 + t44 * t13, t57, -t49 * t7 + t60 * t8, -t2 * r_i_i_C(2) + t62 * t1, t7; -t50 * pkin(3) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t39 * t22 + t8 * t30 + t43 * t31 + t60 * t7 + (-t36 * pkin(5) - pkin(9)) * t13 (-t39 * t38 + t43 * t59) * pkin(2) + t46 * t11 + t44 * t10, -t56, -t49 * t3 + t60 * t4 (t10 * t36 - t4 * t40) * r_i_i_C(2) + t62 * (-t10 * t40 - t4 * t36) t3; 0, t44 * t19 - t46 * t20 + t33 * t61, t34, -t49 * t15 + t60 * t16 (-t16 * t40 + t19 * t36) * r_i_i_C(2) + t62 * (-t16 * t36 - t19 * t40) t15;];
Ja_transl  = t5;
