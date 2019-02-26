% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP6_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP6_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:27
% EndTime: 2019-02-26 20:46:27
% DurationCPUTime: 0.09s
% Computational Cost: add. (87->25), mult. (98->35), div. (0->0), fcn. (111->7), ass. (0->22)
t18 = pkin(3) + pkin(8) + r_i_i_C(3);
t8 = pkin(9) + qJ(3);
t7 = cos(t8);
t16 = t18 * t7;
t6 = sin(t8);
t24 = t16 + t6 * qJ(4) + cos(pkin(9)) * pkin(2) + pkin(1);
t23 = pkin(4) + pkin(7) + qJ(2);
t10 = sin(qJ(5));
t13 = cos(qJ(1));
t22 = t10 * t13;
t11 = sin(qJ(1));
t21 = t11 * t10;
t12 = cos(qJ(5));
t20 = t11 * t12;
t19 = t12 * t13;
t15 = r_i_i_C(1) * t10 + r_i_i_C(2) * t12 + qJ(4);
t14 = t15 * t7 - t18 * t6;
t4 = -t6 * t21 + t19;
t3 = t6 * t20 + t22;
t2 = t6 * t22 + t20;
t1 = t6 * t19 - t21;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t24 * t11 + t23 * t13, t11, t14 * t13, t13 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t23 * t11 + t24 * t13, -t13, t14 * t11, t11 * t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 0, t15 * t6 + t16, -t7 (-r_i_i_C(1) * t12 + r_i_i_C(2) * t10) * t7, 0;];
Ja_transl  = t5;
