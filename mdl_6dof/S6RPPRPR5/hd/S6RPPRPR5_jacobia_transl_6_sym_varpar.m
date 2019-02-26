% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR5_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR5_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:00
% EndTime: 2019-02-26 20:28:00
% DurationCPUTime: 0.10s
% Computational Cost: add. (78->28), mult. (98->38), div. (0->0), fcn. (113->8), ass. (0->21)
t11 = sin(qJ(4));
t13 = cos(qJ(4));
t22 = r_i_i_C(3) + pkin(8) + qJ(5);
t19 = t22 * t13;
t5 = cos(pkin(9)) * pkin(5) + pkin(4);
t23 = -t11 * t5 - pkin(1) - qJ(3) + t19;
t12 = sin(qJ(1));
t21 = t11 * t12;
t14 = cos(qJ(1));
t20 = t11 * t14;
t18 = -pkin(5) * sin(pkin(9)) - pkin(7) + qJ(2);
t8 = pkin(9) + qJ(6);
t6 = sin(t8);
t7 = cos(t8);
t16 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6 + t5;
t15 = t22 * t11 + t16 * t13;
t4 = -t12 * t6 + t7 * t20;
t3 = -t12 * t7 - t6 * t20;
t2 = -t14 * t6 - t7 * t21;
t1 = -t14 * t7 + t6 * t21;
t9 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t12 + t18 * t14, t12, t14, t15 * t14, -t14 * t13, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t18 * t12 - t23 * t14, -t14, t12, t15 * t12, -t12 * t13, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 0, 0, -t16 * t11 + t19, t11 (-r_i_i_C(1) * t6 - r_i_i_C(2) * t7) * t13;];
Ja_transl  = t9;
