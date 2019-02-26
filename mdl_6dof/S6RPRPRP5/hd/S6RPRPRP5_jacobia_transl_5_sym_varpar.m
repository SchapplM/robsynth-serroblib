% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRP5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:52
% EndTime: 2019-02-26 20:45:52
% DurationCPUTime: 0.10s
% Computational Cost: add. (111->28), mult. (96->39), div. (0->0), fcn. (109->9), ass. (0->22)
t12 = pkin(9) + qJ(3);
t10 = cos(t12);
t25 = r_i_i_C(3) + pkin(8) + qJ(4);
t8 = sin(t12);
t22 = t25 * t8;
t5 = cos(pkin(10)) * pkin(4) + pkin(3);
t26 = t22 + t10 * t5 + cos(pkin(9)) * pkin(2) + pkin(1);
t16 = sin(qJ(1));
t24 = t10 * t16;
t17 = cos(qJ(1));
t23 = t10 * t17;
t20 = pkin(4) * sin(pkin(10)) + pkin(7) + qJ(2);
t11 = pkin(10) + qJ(5);
t7 = sin(t11);
t9 = cos(t11);
t19 = r_i_i_C(1) * t9 - r_i_i_C(2) * t7 + t5;
t18 = t25 * t10 - t19 * t8;
t4 = t16 * t7 + t9 * t23;
t3 = t16 * t9 - t7 * t23;
t2 = t17 * t7 - t9 * t24;
t1 = t17 * t9 + t7 * t24;
t6 = [t2 * r_i_i_C(1) + t1 * r_i_i_C(2) - t26 * t16 + t20 * t17, t16, t18 * t17, t17 * t8, r_i_i_C(1) * t3 - r_i_i_C(2) * t4, 0; t4 * r_i_i_C(1) + t3 * r_i_i_C(2) + t20 * t16 + t26 * t17, -t17, t18 * t16, t16 * t8, -r_i_i_C(1) * t1 + r_i_i_C(2) * t2, 0; 0, 0, t19 * t10 + t22, -t10 (-r_i_i_C(1) * t7 - r_i_i_C(2) * t9) * t8, 0;];
Ja_transl  = t6;
