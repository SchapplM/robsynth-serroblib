% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRP2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:30:41
% EndTime: 2019-02-26 20:30:41
% DurationCPUTime: 0.10s
% Computational Cost: add. (111->27), mult. (87->35), div. (0->0), fcn. (97->9), ass. (0->23)
t23 = pkin(8) + r_i_i_C(3);
t10 = pkin(10) + qJ(4);
t6 = sin(t10);
t18 = t23 * t6;
t8 = cos(t10);
t24 = t18 + t8 * pkin(4) + cos(pkin(10)) * pkin(3) + pkin(2);
t14 = cos(qJ(5));
t11 = qJ(1) + pkin(9);
t7 = sin(t11);
t22 = t14 * t7;
t9 = cos(t11);
t21 = t14 * t9;
t13 = sin(qJ(5));
t20 = t7 * t13;
t19 = t9 * t13;
t16 = r_i_i_C(1) * t14 - r_i_i_C(2) * t13 + pkin(4);
t15 = -t16 * t6 + t23 * t8;
t12 = -pkin(7) - qJ(3);
t4 = t8 * t21 + t20;
t3 = -t8 * t19 + t22;
t2 = -t8 * t22 + t19;
t1 = t8 * t20 + t21;
t5 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t9 * t12 - t24 * t7, 0, t7, t15 * t9, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) - t7 * t12 + t24 * t9, 0, -t9, t15 * t7, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, 1, 0, t16 * t8 + t18 (-r_i_i_C(1) * t13 - r_i_i_C(2) * t14) * t6, 0;];
Ja_transl  = t5;
