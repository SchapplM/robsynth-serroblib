% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPPR1_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR1_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobia_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:11
% EndTime: 2019-02-26 20:39:11
% DurationCPUTime: 0.11s
% Computational Cost: add. (160->32), mult. (103->42), div. (0->0), fcn. (116->12), ass. (0->23)
t27 = r_i_i_C(3) + pkin(8) + qJ(5);
t15 = qJ(3) + pkin(10);
t8 = sin(t15);
t30 = cos(qJ(3)) * pkin(3) + t27 * t8;
t11 = cos(t15);
t5 = cos(pkin(11)) * pkin(5) + pkin(4);
t29 = t11 * t5 + pkin(2) + t30;
t16 = qJ(1) + pkin(9);
t9 = sin(t16);
t28 = t11 * t9;
t12 = cos(t16);
t26 = t11 * t12;
t23 = pkin(5) * sin(pkin(11)) + qJ(4) + pkin(7);
t14 = pkin(11) + qJ(6);
t10 = cos(t14);
t7 = sin(t14);
t22 = r_i_i_C(1) * t10 - r_i_i_C(2) * t7 + t5;
t21 = -sin(qJ(3)) * pkin(3) + t27 * t11 - t22 * t8;
t4 = t10 * t26 + t7 * t9;
t3 = t10 * t9 - t7 * t26;
t2 = -t10 * t28 + t12 * t7;
t1 = t10 * t12 + t7 * t28;
t6 = [-sin(qJ(1)) * pkin(1) + t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t12 - t29 * t9, 0, t21 * t12, t9, t12 * t8, r_i_i_C(1) * t3 - r_i_i_C(2) * t4; cos(qJ(1)) * pkin(1) + t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t23 * t9 + t29 * t12, 0, t21 * t9, -t12, t9 * t8, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2; 0, 1, t22 * t11 + t30, 0, -t11 (-r_i_i_C(1) * t7 - r_i_i_C(2) * t10) * t8;];
Ja_transl  = t6;
