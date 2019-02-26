% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:10
% EndTime: 2019-02-26 21:30:11
% DurationCPUTime: 0.21s
% Computational Cost: add. (315->60), mult. (792->92), div. (0->0), fcn. (1045->12), ass. (0->41)
t31 = sin(pkin(11));
t36 = sin(qJ(2));
t40 = cos(qJ(2));
t50 = cos(pkin(11));
t25 = -t40 * t31 - t36 * t50;
t43 = -t36 * t31 + t40 * t50;
t37 = sin(qJ(1));
t41 = cos(qJ(1));
t33 = cos(pkin(6));
t51 = t25 * t33;
t16 = t37 * t51 + t41 * t43;
t63 = -t37 * t43 + t41 * t51;
t35 = sin(qJ(5));
t39 = cos(qJ(5));
t34 = sin(qJ(6));
t38 = cos(qJ(6));
t46 = t38 * r_i_i_C(1) - t34 * r_i_i_C(2) + pkin(5);
t60 = pkin(10) + r_i_i_C(3);
t62 = t46 * t35 - t60 * t39 + qJ(4);
t61 = pkin(3) + pkin(9);
t59 = t40 * pkin(2);
t58 = t33 * t40;
t32 = sin(pkin(6));
t55 = t37 * t32;
t52 = t41 * t32;
t49 = -t33 * t36 * pkin(2) + (pkin(4) + pkin(8) + qJ(3)) * t32;
t22 = t43 * t33;
t12 = t41 * t22 + t37 * t25;
t8 = t12 * t35 + t39 * t52;
t45 = -t12 * t39 + t35 * t52;
t44 = t34 * r_i_i_C(1) + t38 * r_i_i_C(2) + t61;
t30 = pkin(1) + t59;
t21 = t25 * t32;
t20 = t43 * t32;
t18 = -t20 * t35 + t33 * t39;
t15 = -t37 * t22 + t41 * t25;
t4 = -t15 * t35 + t39 * t55;
t3 = t15 * t39 + t35 * t55;
t2 = t16 * t34 + t4 * t38;
t1 = t16 * t38 - t4 * t34;
t5 = [t12 * qJ(4) - t37 * t30 + t49 * t41 + t44 * t63 + t60 * t45 + t46 * t8 (-t41 * t36 - t37 * t58) * pkin(2) + t44 * t15 + t62 * t16, t55, -t15, -t46 * t3 + t60 * t4, t1 * r_i_i_C(1) - t2 * r_i_i_C(2); t4 * pkin(5) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t15 * qJ(4) + t61 * t16 + t60 * t3 + t41 * t30 + t49 * t37 (-t37 * t36 + t41 * t58) * pkin(2) + t44 * t12 - t62 * t63, -t52, -t12, t46 * t45 - t60 * t8 (t8 * t34 - t38 * t63) * r_i_i_C(1) + (t34 * t63 + t8 * t38) * r_i_i_C(2); 0, t44 * t20 - t62 * t21 + t32 * t59, t33, -t20, t60 * t18 + t46 * (-t20 * t39 - t33 * t35) (-t18 * t34 - t21 * t38) * r_i_i_C(1) + (-t18 * t38 + t21 * t34) * r_i_i_C(2);];
Ja_transl  = t5;
