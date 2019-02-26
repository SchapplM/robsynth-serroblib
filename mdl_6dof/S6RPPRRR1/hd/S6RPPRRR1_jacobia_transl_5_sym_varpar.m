% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:34:43
% EndTime: 2019-02-26 20:34:43
% DurationCPUTime: 0.07s
% Computational Cost: add. (87->15), mult. (45->17), div. (0->0), fcn. (47->9), ass. (0->14)
t11 = pkin(11) + qJ(4);
t9 = qJ(5) + t11;
t4 = sin(t9);
t5 = cos(t9);
t16 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
t19 = t16 + pkin(4) * cos(t11);
t17 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(3);
t15 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
t14 = cos(pkin(11)) * pkin(3) + pkin(2) + t19;
t13 = -pkin(4) * sin(t11) + t15;
t12 = qJ(1) + pkin(10);
t8 = cos(t12);
t7 = sin(t12);
t1 = [-sin(qJ(1)) * pkin(1) + t17 * t8 - t14 * t7, 0, t7, t13 * t8, t15 * t8, 0; cos(qJ(1)) * pkin(1) + t17 * t7 + t14 * t8, 0, -t8, t13 * t7, t15 * t7, 0; 0, 1, 0, t19, t16, 0;];
Ja_transl  = t1;
