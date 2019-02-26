% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:49:02
% EndTime: 2019-02-26 20:49:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (90->16), mult. (48->18), div. (0->0), fcn. (50->10), ass. (0->14)
t12 = qJ(3) + pkin(11);
t9 = qJ(5) + t12;
t5 = sin(t9);
t6 = cos(t9);
t17 = t6 * r_i_i_C(1) - t5 * r_i_i_C(2);
t21 = t17 + cos(qJ(3)) * pkin(3) + pkin(4) * cos(t12);
t19 = r_i_i_C(3) + pkin(8) + qJ(4) + pkin(7);
t16 = -r_i_i_C(1) * t5 - r_i_i_C(2) * t6;
t15 = pkin(2) + t21;
t14 = t16 - pkin(4) * sin(t12) - sin(qJ(3)) * pkin(3);
t13 = qJ(1) + pkin(10);
t8 = cos(t13);
t7 = sin(t13);
t1 = [-sin(qJ(1)) * pkin(1) + t19 * t8 - t15 * t7, 0, t14 * t8, t7, t16 * t8, 0; cos(qJ(1)) * pkin(1) + t19 * t7 + t15 * t8, 0, t14 * t7, -t8, t16 * t7, 0; 0, 1, t21, 0, t17, 0;];
Ja_transl  = t1;
