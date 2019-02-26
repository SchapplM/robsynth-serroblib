% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPPR2_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR2_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:22:36
% EndTime: 2019-02-26 21:22:36
% DurationCPUTime: 0.08s
% Computational Cost: add. (80->17), mult. (85->20), div. (0->0), fcn. (97->8), ass. (0->14)
t6 = sin(pkin(10));
t7 = cos(pkin(10));
t15 = r_i_i_C(1) * t6 + r_i_i_C(2) * t7 + qJ(4);
t16 = pkin(3) + r_i_i_C(3) + qJ(5);
t5 = qJ(2) + pkin(9);
t2 = sin(t5);
t3 = cos(t5);
t18 = cos(qJ(2)) * pkin(2) + t15 * t2 + t16 * t3;
t17 = pkin(1) + t18;
t14 = r_i_i_C(1) * t7 - r_i_i_C(2) * t6 + pkin(4) + pkin(7) + qJ(3);
t12 = -sin(qJ(2)) * pkin(2) + t15 * t3 - t16 * t2;
t11 = cos(qJ(1));
t10 = sin(qJ(1));
t1 = [-t17 * t10 + t14 * t11, t12 * t11, t10, t11 * t2, t11 * t3, 0; t14 * t10 + t17 * t11, t12 * t10, -t11, t10 * t2, t10 * t3, 0; 0, t18, 0, -t3, t2, 0;];
Ja_transl  = t1;
