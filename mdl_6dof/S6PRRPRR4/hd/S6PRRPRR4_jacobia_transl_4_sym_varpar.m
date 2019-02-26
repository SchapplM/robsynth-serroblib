% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function Ja_transl = S6PRRPRR4_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR4_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:52
% EndTime: 2019-02-26 20:05:52
% DurationCPUTime: 0.11s
% Computational Cost: add. (66->21), mult. (169->39), div. (0->0), fcn. (211->8), ass. (0->22)
t28 = pkin(3) + r_i_i_C(1);
t27 = pkin(8) + r_i_i_C(2);
t14 = sin(pkin(6));
t17 = sin(qJ(3));
t26 = t14 * t17;
t19 = cos(qJ(3));
t25 = t14 * t19;
t16 = cos(pkin(6));
t18 = sin(qJ(2));
t24 = t16 * t18;
t20 = cos(qJ(2));
t23 = t16 * t20;
t22 = r_i_i_C(3) + qJ(4);
t21 = t22 * t17 + t19 * t28 + pkin(2);
t15 = cos(pkin(11));
t13 = sin(pkin(11));
t9 = -t16 * t19 + t18 * t26;
t8 = -t13 * t24 + t15 * t20;
t6 = t13 * t20 + t15 * t24;
t3 = -t13 * t25 + t17 * t8;
t1 = t15 * t25 + t17 * t6;
t2 = [0, t27 * t8 + t21 * (-t13 * t23 - t15 * t18) t22 * (t13 * t26 + t19 * t8) - t28 * t3, t3, 0, 0; 0, t27 * t6 + t21 * (-t13 * t18 + t15 * t23) t22 * (-t15 * t26 + t19 * t6) - t28 * t1, t1, 0, 0; 1 (t18 * t27 + t20 * t21) * t14, -t28 * t9 + t22 * (t16 * t17 + t18 * t25) t9, 0, 0;];
Ja_transl  = t2;
