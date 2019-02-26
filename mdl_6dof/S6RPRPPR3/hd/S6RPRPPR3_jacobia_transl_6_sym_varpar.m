% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR3_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR3_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:40:24
% EndTime: 2019-02-26 20:40:24
% DurationCPUTime: 0.09s
% Computational Cost: add. (110->27), mult. (116->38), div. (0->0), fcn. (129->8), ass. (0->21)
t11 = cos(qJ(3));
t15 = pkin(3) + pkin(4) + pkin(8) + r_i_i_C(3);
t14 = t15 * t11;
t17 = pkin(5) + qJ(4);
t9 = sin(qJ(3));
t21 = -t17 * t9 - pkin(2) - t14;
t8 = sin(qJ(6));
t19 = t8 * t9;
t10 = cos(qJ(6));
t18 = t10 * t9;
t16 = pkin(7) - qJ(5);
t13 = r_i_i_C(1) * t10 - r_i_i_C(2) * t8 + t17;
t12 = t13 * t11 - t15 * t9;
t7 = qJ(1) + pkin(9);
t6 = cos(t7);
t5 = sin(t7);
t4 = t6 * t18 - t5 * t8;
t3 = -t10 * t5 - t6 * t19;
t2 = -t5 * t18 - t6 * t8;
t1 = -t10 * t6 + t5 * t19;
t20 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t16 * t6 + t21 * t5, 0, t12 * t6, t6 * t9, -t5, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t16 * t5 - t21 * t6, 0, t12 * t5, t5 * t9, t6, -t1 * r_i_i_C(1) + r_i_i_C(2) * t2; 0, 1, t13 * t9 + t14, -t11, 0 (r_i_i_C(1) * t8 + r_i_i_C(2) * t10) * t11;];
Ja_transl  = t20;
