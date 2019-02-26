% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR7_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR7_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:42:32
% EndTime: 2019-02-26 20:42:33
% DurationCPUTime: 0.09s
% Computational Cost: add. (92->27), mult. (107->36), div. (0->0), fcn. (122->8), ass. (0->22)
t17 = pkin(4) + pkin(8) + r_i_i_C(3);
t7 = qJ(3) + pkin(9);
t5 = sin(t7);
t26 = -t17 * t5 - sin(qJ(3)) * pkin(3);
t12 = cos(qJ(6));
t9 = sin(qJ(6));
t16 = r_i_i_C(1) * t9 + r_i_i_C(2) * t12 + qJ(5);
t6 = cos(t7);
t25 = t16 * t5 + t17 * t6 + cos(qJ(3)) * pkin(3);
t11 = sin(qJ(1));
t22 = t11 * t9;
t14 = cos(qJ(1));
t21 = t14 * t9;
t20 = t11 * t12;
t19 = t12 * t14;
t18 = pkin(1) + pkin(5) + qJ(4) + pkin(7);
t15 = -t6 * qJ(5) + qJ(2) - t26;
t4 = -t6 * t22 + t19;
t3 = -t6 * t20 - t21;
t2 = -t6 * t21 - t20;
t1 = -t6 * t19 + t22;
t8 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t18 * t11 + t15 * t14, t11, t25 * t11, t14, -t11 * t6, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t15 * t11 + t18 * t14, -t14, -t25 * t14, t11, t14 * t6, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 0, t16 * t6 + t26, 0, t5 (r_i_i_C(1) * t12 - r_i_i_C(2) * t9) * t5;];
Ja_transl  = t8;
