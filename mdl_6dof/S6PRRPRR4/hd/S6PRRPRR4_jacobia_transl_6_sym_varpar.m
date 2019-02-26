% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPRR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:52
% EndTime: 2019-02-26 20:05:52
% DurationCPUTime: 0.20s
% Computational Cost: add. (309->51), mult. (815->88), div. (0->0), fcn. (1061->12), ass. (0->42)
t35 = sin(qJ(5));
t36 = sin(qJ(3));
t39 = cos(qJ(5));
t40 = cos(qJ(3));
t34 = sin(qJ(6));
t38 = cos(qJ(6));
t46 = t38 * r_i_i_C(1) - t34 * r_i_i_C(2) + pkin(5);
t61 = r_i_i_C(3) + pkin(10);
t67 = t61 * (t35 * t40 - t36 * t39) + (t35 * t36 + t39 * t40) * t46;
t62 = pkin(3) + pkin(4);
t44 = t36 * qJ(4) + t62 * t40 + pkin(2);
t66 = t44 + t67;
t37 = sin(qJ(2));
t57 = cos(pkin(6));
t32 = sin(pkin(6));
t60 = t32 * t36;
t27 = t37 * t60 - t57 * t40;
t59 = t32 * t40;
t28 = t57 * t36 + t37 * t59;
t16 = t27 * t35 + t28 * t39;
t65 = t46 * (t27 * t39 - t28 * t35) + t61 * t16;
t33 = cos(pkin(11));
t41 = cos(qJ(2));
t56 = sin(pkin(11));
t47 = t57 * t56;
t26 = t33 * t41 - t37 * t47;
t52 = t32 * t56;
t19 = t26 * t36 - t40 * t52;
t20 = t26 * t40 + t36 * t52;
t8 = t19 * t35 + t20 * t39;
t64 = t46 * (t19 * t39 - t20 * t35) + t61 * t8;
t51 = t33 * t57;
t24 = t37 * t51 + t56 * t41;
t17 = t24 * t36 + t33 * t59;
t18 = t24 * t40 - t33 * t60;
t4 = t17 * t35 + t18 * t39;
t63 = t46 * (t17 * t39 - t18 * t35) + t61 * t4;
t58 = t32 * t41;
t45 = -t34 * r_i_i_C(1) - t38 * r_i_i_C(2) + pkin(8) - pkin(9);
t25 = -t33 * t37 - t41 * t47;
t23 = -t56 * t37 + t41 * t51;
t1 = [0, t66 * t25 + t45 * t26, t20 * qJ(4) - t62 * t19 - t64, t19, t64 (t25 * t38 - t8 * t34) * r_i_i_C(1) + (-t25 * t34 - t8 * t38) * r_i_i_C(2); 0, t66 * t23 + t45 * t24, t18 * qJ(4) - t62 * t17 - t63, t17, t63 (t23 * t38 - t4 * t34) * r_i_i_C(1) + (-t23 * t34 - t4 * t38) * r_i_i_C(2); 1 (t45 * t37 + t44 * t41) * t32 + t67 * t58, t28 * qJ(4) - t62 * t27 - t65, t27, t65 (-t16 * t34 + t38 * t58) * r_i_i_C(1) + (-t16 * t38 - t34 * t58) * r_i_i_C(2);];
Ja_transl  = t1;
