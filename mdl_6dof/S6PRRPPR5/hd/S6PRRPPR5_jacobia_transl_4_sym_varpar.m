% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRRPPR5_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR5_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:37
% EndTime: 2019-02-26 20:00:37
% DurationCPUTime: 0.11s
% Computational Cost: add. (66->21), mult. (169->39), div. (0->0), fcn. (211->8), ass. (0->22)
t26 = pkin(3) - r_i_i_C(2);
t25 = pkin(8) + r_i_i_C(1);
t12 = sin(pkin(6));
t15 = sin(qJ(3));
t24 = t12 * t15;
t17 = cos(qJ(3));
t23 = t12 * t17;
t14 = cos(pkin(6));
t16 = sin(qJ(2));
t22 = t14 * t16;
t18 = cos(qJ(2));
t21 = t14 * t18;
t20 = r_i_i_C(3) + qJ(4);
t19 = t20 * t15 + t17 * t26 + pkin(2);
t13 = cos(pkin(10));
t11 = sin(pkin(10));
t9 = -t14 * t17 + t16 * t24;
t8 = -t11 * t22 + t13 * t18;
t6 = t11 * t18 + t13 * t22;
t3 = -t11 * t23 + t15 * t8;
t1 = t13 * t23 + t15 * t6;
t2 = [0, t25 * t8 + t19 * (-t11 * t21 - t13 * t16) t20 * (t11 * t24 + t17 * t8) - t26 * t3, t3, 0, 0; 0, t25 * t6 + t19 * (-t11 * t16 + t13 * t21) t20 * (-t13 * t24 + t17 * t6) - t26 * t1, t1, 0, 0; 1 (t16 * t25 + t18 * t19) * t12, -t26 * t9 + t20 * (t14 * t15 + t16 * t23) t9, 0, 0;];
Ja_transl  = t2;
