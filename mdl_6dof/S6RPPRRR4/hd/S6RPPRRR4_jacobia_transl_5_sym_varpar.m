% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR4_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR4_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:36:28
% EndTime: 2019-02-26 20:36:28
% DurationCPUTime: 0.10s
% Computational Cost: add. (78->28), mult. (161->40), div. (0->0), fcn. (207->8), ass. (0->22)
t11 = cos(qJ(4));
t10 = cos(qJ(5));
t8 = sin(qJ(5));
t13 = r_i_i_C(1) * t10 - r_i_i_C(2) * t8 + pkin(4);
t22 = pkin(8) + r_i_i_C(3);
t9 = sin(qJ(4));
t15 = t22 * t9;
t25 = -t13 * t11 - t15;
t24 = pkin(1) + pkin(2);
t21 = t11 * t8;
t20 = cos(qJ(1));
t19 = sin(qJ(1));
t18 = t10 * t11;
t17 = cos(pkin(10));
t16 = sin(pkin(10));
t14 = t8 * r_i_i_C(1) + t10 * r_i_i_C(2);
t12 = -t22 * t11 + t13 * t9;
t4 = t20 * t16 - t19 * t17;
t3 = -t19 * t16 - t20 * t17;
t2 = -t3 * t18 + t4 * t8;
t1 = t4 * t10 + t3 * t21;
t5 = [t20 * qJ(2) + (pkin(7) + t14) * t3 + (pkin(3) - t25) * t4 - t24 * t19, t19, 0, t12 * t3, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t19 * qJ(2) + t4 * pkin(7) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + (-t11 * pkin(4) - pkin(3) - t15) * t3 + t24 * t20, -t20, 0, t12 * t4 (-t3 * t10 + t4 * t21) * r_i_i_C(1) + (t4 * t18 + t3 * t8) * r_i_i_C(2), 0; 0, 0, -1, t25, t14 * t9, 0;];
Ja_transl  = t5;
