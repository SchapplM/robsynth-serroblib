% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobia_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:32:34
% EndTime: 2019-02-26 21:32:34
% DurationCPUTime: 0.18s
% Computational Cost: add. (176->43), mult. (315->66), div. (0->0), fcn. (381->10), ass. (0->33)
t24 = cos(qJ(2));
t21 = sin(qJ(2));
t33 = r_i_i_C(3) - qJ(3) + pkin(9) + pkin(8);
t30 = t33 * t21;
t42 = -t24 * pkin(2) - pkin(1) + t30;
t18 = sin(pkin(10));
t19 = cos(pkin(10));
t17 = qJ(5) + qJ(6);
t15 = sin(t17);
t16 = cos(t17);
t20 = sin(qJ(5));
t31 = t20 * pkin(5) + qJ(4);
t28 = t15 * r_i_i_C(1) + t16 * r_i_i_C(2) + t31;
t23 = cos(qJ(5));
t37 = t23 * pkin(5) + pkin(3) + pkin(4);
t29 = t16 * r_i_i_C(1) - t15 * r_i_i_C(2) + t37;
t41 = t28 * t18 + t29 * t19 + pkin(2);
t25 = cos(qJ(1));
t35 = t25 * t18;
t22 = sin(qJ(1));
t36 = t22 * t24;
t10 = t19 * t36 - t35;
t34 = t25 * t19;
t9 = t18 * t36 + t34;
t40 = (-t10 * t15 + t9 * t16) * r_i_i_C(1) + (-t10 * t16 - t9 * t15) * r_i_i_C(2);
t11 = -t22 * t19 + t24 * t35;
t12 = t22 * t18 + t24 * t34;
t5 = t11 * t16 - t12 * t15;
t6 = t11 * t15 + t12 * t16;
t39 = t5 * r_i_i_C(1) - t6 * r_i_i_C(2);
t38 = ((-t15 * t19 + t16 * t18) * r_i_i_C(1) + (-t15 * t18 - t16 * t19) * r_i_i_C(2)) * t21;
t27 = -t41 * t21 - t33 * t24;
t1 = [t25 * pkin(7) - t29 * t10 + t42 * t22 - t28 * t9, t27 * t25, t25 * t21, t11 (t11 * t23 - t12 * t20) * pkin(5) + t39, t39; t22 * pkin(7) + t6 * r_i_i_C(1) + t5 * r_i_i_C(2) + t31 * t11 + t37 * t12 - t42 * t25, t27 * t22, t22 * t21, t9 (-t10 * t20 + t23 * t9) * pkin(5) + t40, t40; 0, t41 * t24 - t30, -t24, t21 * t18 (t18 * t23 - t19 * t20) * t21 * pkin(5) + t38, t38;];
Ja_transl  = t1;
