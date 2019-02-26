% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPRPR7_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR7_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:09
% EndTime: 2019-02-26 20:29:09
% DurationCPUTime: 0.08s
% Computational Cost: add. (64->17), mult. (71->17), div. (0->0), fcn. (82->7), ass. (0->14)
t3 = pkin(9) + qJ(4);
t1 = sin(t3);
t4 = sin(pkin(10));
t6 = cos(pkin(10));
t12 = r_i_i_C(1) * t6 - r_i_i_C(2) * t4 + pkin(4);
t13 = r_i_i_C(3) + qJ(5);
t2 = cos(t3);
t15 = t13 * t1 + t12 * t2;
t14 = -t12 * t1 + t13 * t2;
t11 = r_i_i_C(1) * t4 + r_i_i_C(2) * t6 + pkin(1) + pkin(7) + qJ(3);
t10 = sin(pkin(9)) * pkin(3) + qJ(2) - t14;
t9 = cos(qJ(1));
t8 = sin(qJ(1));
t5 = [t10 * t9 - t11 * t8, t8, t9, t15 * t8, -t8 * t2, 0; t10 * t8 + t11 * t9, -t9, t8, -t15 * t9, t9 * t2, 0; 0, 0, 0, t14, t1, 0;];
Ja_transl  = t5;
