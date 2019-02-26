% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRRRPR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:30:46
% EndTime: 2019-02-26 22:30:46
% DurationCPUTime: 0.08s
% Computational Cost: add. (144->16), mult. (74->20), div. (0->0), fcn. (76->10), ass. (0->16)
t15 = qJ(2) + qJ(3);
t13 = qJ(4) + t15;
t9 = pkin(11) + t13;
t5 = sin(t9);
t6 = cos(t9);
t24 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2) + pkin(4) * cos(t13);
t22 = t24 + pkin(3) * cos(t15);
t28 = cos(qJ(2)) * pkin(2) + t22;
t18 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6 - pkin(4) * sin(t13);
t19 = -pkin(3) * sin(t15) + t18;
t25 = r_i_i_C(3) + qJ(5) + pkin(9) + pkin(8) + pkin(7);
t21 = pkin(1) + t28;
t20 = -sin(qJ(2)) * pkin(2) + t19;
t17 = cos(qJ(1));
t16 = sin(qJ(1));
t1 = [-t21 * t16 + t25 * t17, t20 * t17, t19 * t17, t18 * t17, t16, 0; t25 * t16 + t21 * t17, t20 * t16, t19 * t16, t18 * t16, -t17, 0; 0, t28, t22, t24, 0, 0;];
Ja_transl  = t1;
