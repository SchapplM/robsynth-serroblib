% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_jacobia_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:43:07
% EndTime: 2019-02-26 20:43:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (58->27), mult. (116->34), div. (0->0), fcn. (131->6), ass. (0->20)
t13 = pkin(3) + pkin(4) + pkin(8) + r_i_i_C(3);
t6 = sin(qJ(3));
t19 = t13 * t6;
t15 = -pkin(5) - qJ(4);
t5 = sin(qJ(6));
t8 = cos(qJ(6));
t12 = r_i_i_C(1) * t8 - r_i_i_C(2) * t5 - t15;
t9 = cos(qJ(3));
t18 = t12 * t6 + t13 * t9;
t7 = sin(qJ(1));
t17 = t7 * t9;
t10 = cos(qJ(1));
t16 = t10 * t9;
t14 = pkin(1) + pkin(7) - qJ(5);
t11 = t15 * t9 + qJ(2) + t19;
t4 = -t8 * t16 + t7 * t5;
t3 = t5 * t16 + t7 * t8;
t2 = t10 * t5 + t8 * t17;
t1 = -t10 * t8 + t5 * t17;
t20 = [t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t11 * t10 - t14 * t7, t7, t18 * t7, -t17, -t10, r_i_i_C(1) * t1 + r_i_i_C(2) * t2; -t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t14 * t10 + t11 * t7, -t10, -t18 * t10, t16, -t7, -r_i_i_C(1) * t3 + r_i_i_C(2) * t4; 0, 0, t12 * t9 - t19, t6, 0 (-r_i_i_C(1) * t5 - r_i_i_C(2) * t8) * t6;];
Ja_transl  = t20;
