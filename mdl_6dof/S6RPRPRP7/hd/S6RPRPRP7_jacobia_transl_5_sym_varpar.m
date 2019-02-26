% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP7
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

function Ja_transl = S6RPRPRP7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:56
% EndTime: 2019-02-26 20:46:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (76->26), mult. (94->34), div. (0->0), fcn. (106->8), ass. (0->22)
t23 = pkin(8) + r_i_i_C(3);
t7 = qJ(3) + pkin(9);
t6 = cos(t7);
t26 = t23 * t6 - sin(qJ(3)) * pkin(3);
t12 = cos(qJ(5));
t9 = sin(qJ(5));
t16 = r_i_i_C(1) * t12 - r_i_i_C(2) * t9 + pkin(4);
t5 = sin(t7);
t25 = t16 * t6 + t23 * t5 + cos(qJ(3)) * pkin(3);
t24 = pkin(1) + qJ(4) + pkin(7);
t11 = sin(qJ(1));
t20 = t11 * t9;
t14 = cos(qJ(1));
t19 = t14 * t9;
t18 = t11 * t12;
t17 = t12 * t14;
t15 = t5 * pkin(4) + qJ(2) - t26;
t4 = t5 * t17 - t20;
t3 = t5 * t19 + t18;
t2 = t5 * t18 + t19;
t1 = -t5 * t20 + t17;
t8 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t24 * t11 + t15 * t14, t11, t25 * t11, t14, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t11 + t24 * t14, -t14, -t25 * t14, t11, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, 0, -t16 * t5 + t26, 0 (-r_i_i_C(1) * t9 - r_i_i_C(2) * t12) * t6, 0;];
Ja_transl  = t8;
