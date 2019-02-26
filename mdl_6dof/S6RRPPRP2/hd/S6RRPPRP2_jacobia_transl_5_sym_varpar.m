% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRP2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:25:21
% EndTime: 2019-02-26 21:25:21
% DurationCPUTime: 0.10s
% Computational Cost: add. (90->26), mult. (103->36), div. (0->0), fcn. (116->8), ass. (0->22)
t20 = pkin(3) + pkin(8) + r_i_i_C(3);
t9 = qJ(2) + pkin(9);
t7 = cos(t9);
t27 = t20 * t7 + cos(qJ(2)) * pkin(2);
t6 = sin(t9);
t26 = t6 * qJ(4) + pkin(1) + t27;
t25 = pkin(4) + qJ(3) + pkin(7);
t11 = sin(qJ(5));
t15 = cos(qJ(1));
t24 = t11 * t15;
t13 = sin(qJ(1));
t23 = t13 * t11;
t14 = cos(qJ(5));
t22 = t13 * t14;
t21 = t14 * t15;
t17 = r_i_i_C(1) * t11 + r_i_i_C(2) * t14 + qJ(4);
t16 = -sin(qJ(2)) * pkin(2) + t17 * t7 - t20 * t6;
t4 = -t6 * t23 + t21;
t3 = t6 * t22 + t24;
t2 = t6 * t24 + t22;
t1 = t6 * t21 - t23;
t5 = [t4 * r_i_i_C(1) - t3 * r_i_i_C(2) - t26 * t13 + t25 * t15, t16 * t15, t13, t15 * t6, r_i_i_C(1) * t1 - r_i_i_C(2) * t2, 0; t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t25 * t13 + t26 * t15, t16 * t13, -t15, t13 * t6, r_i_i_C(1) * t3 + r_i_i_C(2) * t4, 0; 0, t17 * t6 + t27, 0, -t7 (-r_i_i_C(1) * t14 + r_i_i_C(2) * t11) * t7, 0;];
Ja_transl  = t5;
