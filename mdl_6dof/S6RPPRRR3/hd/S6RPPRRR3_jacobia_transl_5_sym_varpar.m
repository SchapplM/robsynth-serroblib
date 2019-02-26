% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:35:56
% EndTime: 2019-02-26 20:35:57
% DurationCPUTime: 0.09s
% Computational Cost: add. (80->26), mult. (87->36), div. (0->0), fcn. (97->8), ass. (0->20)
t11 = cos(qJ(4));
t15 = pkin(8) + r_i_i_C(3);
t19 = t15 * t11;
t10 = cos(qJ(5));
t8 = sin(qJ(5));
t13 = r_i_i_C(1) * t10 - r_i_i_C(2) * t8 + pkin(4);
t9 = sin(qJ(4));
t18 = t13 * t11 + t15 * t9;
t17 = t8 * t9;
t16 = pkin(2) + pkin(7);
t14 = t10 * t9;
t12 = t9 * pkin(4) + qJ(3) - t19;
t7 = qJ(1) + pkin(10);
t6 = cos(t7);
t5 = sin(t7);
t4 = t6 * t14 - t5 * t8;
t3 = t10 * t5 + t6 * t17;
t2 = t5 * t14 + t6 * t8;
t1 = t10 * t6 - t5 * t17;
t20 = [-sin(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t16 * t5 + t12 * t6, 0, t5, t18 * t5, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; cos(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t6 + t12 * t5, 0, -t6, -t18 * t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 1, 0, -t13 * t9 + t19 (-r_i_i_C(1) * t8 - r_i_i_C(2) * t10) * t11, 0;];
Ja_transl  = t20;
