% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRRPR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:50
% EndTime: 2019-02-26 21:00:50
% DurationCPUTime: 0.08s
% Computational Cost: add. (98->16), mult. (53->18), div. (0->0), fcn. (55->10), ass. (0->14)
t14 = qJ(3) + qJ(4);
t9 = pkin(11) + t14;
t4 = sin(t9);
t5 = cos(t9);
t19 = t5 * r_i_i_C(1) - t4 * r_i_i_C(2) + pkin(4) * cos(t14);
t23 = cos(qJ(3)) * pkin(3) + t19;
t15 = -r_i_i_C(1) * t4 - r_i_i_C(2) * t5 - pkin(4) * sin(t14);
t20 = r_i_i_C(3) + qJ(5) + pkin(8) + pkin(7);
t17 = pkin(2) + t23;
t16 = -sin(qJ(3)) * pkin(3) + t15;
t13 = qJ(1) + pkin(10);
t8 = cos(t13);
t7 = sin(t13);
t1 = [-sin(qJ(1)) * pkin(1) + t20 * t8 - t17 * t7, 0, t16 * t8, t15 * t8, t7, 0; cos(qJ(1)) * pkin(1) + t20 * t7 + t17 * t8, 0, t16 * t7, t15 * t7, -t8, 0; 0, 1, t23, t19, 0, 0;];
Ja_transl  = t1;
