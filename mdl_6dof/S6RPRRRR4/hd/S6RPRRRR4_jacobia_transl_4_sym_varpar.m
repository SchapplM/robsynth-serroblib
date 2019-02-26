% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRRR4_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR4_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:25
% EndTime: 2019-02-26 21:16:26
% DurationCPUTime: 0.07s
% Computational Cost: add. (63->12), mult. (43->15), div. (0->0), fcn. (45->7), ass. (0->13)
t9 = pkin(11) + qJ(3);
t7 = qJ(4) + t9;
t4 = sin(t7);
t5 = cos(t7);
t15 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2);
t18 = t15 + pkin(3) * cos(t9);
t16 = r_i_i_C(3) + pkin(8) + pkin(7) + qJ(2);
t14 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5;
t13 = cos(pkin(11)) * pkin(2) + pkin(1) + t18;
t12 = -pkin(3) * sin(t9) + t14;
t11 = cos(qJ(1));
t10 = sin(qJ(1));
t1 = [-t10 * t13 + t11 * t16, t10, t12 * t11, t14 * t11, 0, 0; t10 * t16 + t11 * t13, -t11, t12 * t10, t14 * t10, 0, 0; 0, 0, t18, t15, 0, 0;];
Ja_transl  = t1;
