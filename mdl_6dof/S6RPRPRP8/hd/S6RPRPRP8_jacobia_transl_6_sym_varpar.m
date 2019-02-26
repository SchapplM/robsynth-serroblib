% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP8_jacobia_transl_6_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_jacobia_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP8_jacobia_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_jacobia_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:47:36
% EndTime: 2019-02-26 20:47:36
% DurationCPUTime: 0.10s
% Computational Cost: add. (121->29), mult. (153->35), div. (0->0), fcn. (178->8), ass. (0->24)
t25 = pkin(8) + r_i_i_C(2);
t10 = sin(qJ(5));
t13 = cos(qJ(5));
t18 = r_i_i_C(3) + qJ(6);
t26 = pkin(5) + r_i_i_C(1);
t28 = t18 * t10 + t26 * t13 + pkin(4);
t8 = qJ(3) + pkin(9);
t6 = sin(t8);
t7 = cos(t8);
t31 = t28 * t7 + t25 * t6 + cos(qJ(3)) * pkin(3);
t30 = t25 * t7 - sin(qJ(3)) * pkin(3);
t27 = pkin(1) + qJ(4) + pkin(7);
t12 = sin(qJ(1));
t22 = t12 * t10;
t21 = t12 * t13;
t15 = cos(qJ(1));
t20 = t15 * t10;
t19 = t15 * t13;
t16 = t6 * pkin(4) + qJ(2) - t30;
t4 = t6 * t19 - t22;
t3 = t6 * t20 + t21;
t2 = t6 * t21 + t20;
t1 = t6 * t22 - t19;
t5 = [-t27 * t12 + t16 * t15 + t18 * t3 + t26 * t4, t12, t31 * t12, t15, -t26 * t1 + t18 * t2, t1; t18 * t1 + t16 * t12 + t27 * t15 + t26 * t2, -t15, -t31 * t15, t12, -t18 * t4 + t26 * t3, -t3; 0, 0, -t28 * t6 + t30, 0 (-t26 * t10 + t18 * t13) * t7, t7 * t10;];
Ja_transl  = t5;
