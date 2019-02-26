% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPPRR4_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR4_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:24:17
% EndTime: 2019-02-26 20:24:17
% DurationCPUTime: 0.10s
% Computational Cost: add. (84->31), mult. (171->40), div. (0->0), fcn. (223->8), ass. (0->23)
t25 = pkin(1) + pkin(2);
t11 = cos(qJ(5));
t10 = cos(qJ(6));
t8 = sin(qJ(6));
t13 = r_i_i_C(1) * t10 - r_i_i_C(2) * t8 + pkin(5);
t21 = pkin(8) + r_i_i_C(3);
t9 = sin(qJ(5));
t24 = t13 * t11 + t21 * t9;
t23 = t8 * t9;
t22 = -pkin(3) - pkin(7);
t20 = t10 * t9;
t19 = cos(qJ(1));
t18 = sin(qJ(1));
t17 = cos(pkin(9));
t16 = sin(pkin(9));
t15 = t21 * t11;
t14 = t8 * r_i_i_C(1) + t10 * r_i_i_C(2);
t12 = t13 * t9 - t15;
t4 = t19 * t16 - t18 * t17;
t3 = -t18 * t16 - t19 * t17;
t2 = t4 * t20 - t3 * t8;
t1 = -t10 * t3 - t4 * t23;
t5 = [t19 * qJ(2) + (t14 - t22) * t4 + (qJ(4) + t12) * t3 - t25 * t18, t18, 0, t4, t24 * t4, r_i_i_C(1) * t1 - r_i_i_C(2) * t2; t18 * qJ(2) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t22 * t3 + (t9 * pkin(5) + qJ(4) - t15) * t4 + t25 * t19, -t19, 0, -t3, -t24 * t3 (-t10 * t4 + t3 * t23) * r_i_i_C(1) + (t3 * t20 + t4 * t8) * r_i_i_C(2); 0, 0, -1, 0, t12, t14 * t11;];
Ja_transl  = t5;
