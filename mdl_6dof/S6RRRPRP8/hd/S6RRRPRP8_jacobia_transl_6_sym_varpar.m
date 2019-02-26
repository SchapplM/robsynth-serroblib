% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP8
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRPRP8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:13:08
% EndTime: 2019-02-26 22:13:09
% DurationCPUTime: 0.20s
% Computational Cost: add. (156->45), mult. (361->67), div. (0->0), fcn. (435->8), ass. (0->33)
t15 = sin(qJ(3));
t19 = cos(qJ(3));
t18 = cos(qJ(5));
t14 = sin(qJ(5));
t37 = t18 * pkin(5) + pkin(3) + pkin(4);
t29 = t14 * r_i_i_C(2) - t37;
t24 = t18 * r_i_i_C(1) - t29;
t39 = pkin(5) + r_i_i_C(1);
t26 = -t14 * t39 - qJ(4);
t41 = -t18 * r_i_i_C(2) + t26;
t44 = t41 * t15 - t24 * t19 - pkin(2);
t20 = cos(qJ(2));
t16 = sin(qJ(2));
t33 = -r_i_i_C(3) - qJ(6) - pkin(9) + pkin(8);
t30 = t33 * t16;
t43 = t20 * pkin(2) + pkin(1) + t30;
t42 = (t14 * t19 - t15 * t18) * t16;
t21 = cos(qJ(1));
t34 = t21 * t19;
t17 = sin(qJ(1));
t36 = t17 * t20;
t6 = t15 * t36 + t34;
t3 = t6 * t18;
t35 = t21 * t15;
t31 = t14 * pkin(5) + qJ(4);
t28 = -r_i_i_C(1) * t42 - (t14 * t15 + t18 * t19) * t16 * r_i_i_C(2);
t8 = -t17 * t19 + t20 * t35;
t9 = t17 * t15 + t20 * t34;
t1 = -t9 * t14 + t8 * t18;
t22 = t44 * t16 + t33 * t20;
t7 = t19 * t36 - t35;
t2 = t8 * t14 + t9 * t18;
t4 = [t21 * pkin(7) - t3 * r_i_i_C(2) - t43 * t17 - t24 * t7 + t26 * t6, t22 * t21, -t24 * t8 - t41 * t9, t8, -t2 * r_i_i_C(2) + t39 * t1, -t21 * t16; t17 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t43 * t21 + t31 * t8 + t37 * t9, t22 * t17, -t3 * r_i_i_C(1) + t29 * t6 - t41 * t7, t6 (-t6 * t14 - t7 * t18) * r_i_i_C(2) + t39 * (-t7 * t14 + t3) -t17 * t16; 0, -t44 * t20 + t30 (-t15 * t37 + t19 * t31) * t16 - t28, t16 * t15, -pkin(5) * t42 + t28, t20;];
Ja_transl  = t4;
