% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPRPR8_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR8_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:47
% EndTime: 2019-02-26 21:41:47
% DurationCPUTime: 0.11s
% Computational Cost: add. (127->26), mult. (153->39), div. (0->0), fcn. (177->8), ass. (0->23)
t15 = cos(qJ(2));
t13 = sin(qJ(2));
t24 = r_i_i_C(2) + pkin(8) + qJ(3);
t18 = t24 * t13;
t7 = cos(pkin(10)) * pkin(3) + pkin(2);
t27 = t15 * t7 + pkin(1) + t18;
t21 = r_i_i_C(3) + qJ(5);
t25 = pkin(4) + r_i_i_C(1);
t10 = pkin(10) + qJ(4);
t8 = sin(t10);
t9 = cos(t10);
t26 = t21 * t8 + t25 * t9 + t7;
t14 = sin(qJ(1));
t23 = t14 * t15;
t16 = cos(qJ(1));
t22 = t15 * t16;
t19 = pkin(3) * sin(pkin(10)) + pkin(7);
t17 = -t26 * t13 + t24 * t15;
t4 = t14 * t8 + t9 * t22;
t3 = -t14 * t9 + t8 * t22;
t2 = -t16 * t8 + t9 * t23;
t1 = t16 * t9 + t8 * t23;
t5 = [-t21 * t1 - t27 * t14 + t19 * t16 - t25 * t2, t17 * t16, t16 * t13, t21 * t4 - t25 * t3, t3, 0; t19 * t14 + t27 * t16 + t21 * t3 + t25 * t4, t17 * t14, t14 * t13, -t25 * t1 + t21 * t2, t1, 0; 0, t26 * t15 + t18, -t15 (t21 * t9 - t25 * t8) * t13, t13 * t8, 0;];
Ja_transl  = t5;
