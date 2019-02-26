% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR1_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR1_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobia_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:25:29
% EndTime: 2019-02-26 20:25:29
% DurationCPUTime: 0.09s
% Computational Cost: add. (96->18), mult. (69->19), div. (0->0), fcn. (78->9), ass. (0->15)
t8 = sin(pkin(11));
t9 = cos(pkin(11));
t14 = r_i_i_C(1) * t9 - r_i_i_C(2) * t8 + pkin(4);
t15 = r_i_i_C(3) + qJ(5);
t6 = pkin(10) + qJ(4);
t2 = sin(t6);
t4 = cos(t6);
t12 = t14 * t4 + t15 * t2;
t16 = cos(pkin(10)) * pkin(3) + pkin(2) + t12;
t13 = r_i_i_C(1) * t8 + r_i_i_C(2) * t9 + pkin(7) + qJ(3);
t11 = -t14 * t2 + t15 * t4;
t7 = qJ(1) + pkin(9);
t5 = cos(t7);
t3 = sin(t7);
t1 = [-sin(qJ(1)) * pkin(1) + t13 * t5 - t16 * t3, 0, t3, t11 * t5, t5 * t2, 0; cos(qJ(1)) * pkin(1) + t13 * t3 + t16 * t5, 0, -t5, t11 * t3, t3 * t2, 0; 0, 1, 0, t12, -t4, 0;];
Ja_transl  = t1;
