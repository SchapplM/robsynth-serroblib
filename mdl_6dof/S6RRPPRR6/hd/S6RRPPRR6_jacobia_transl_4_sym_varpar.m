% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RRPPRR6_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR6_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_jacobia_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:25
% EndTime: 2019-02-26 21:31:25
% DurationCPUTime: 0.08s
% Computational Cost: add. (40->24), mult. (83->34), div. (0->0), fcn. (98->6), ass. (0->17)
t16 = pkin(2) + pkin(3);
t7 = sin(qJ(2));
t9 = cos(qJ(2));
t17 = t7 * qJ(3) + t16 * t9 + pkin(1);
t15 = pkin(7) - r_i_i_C(3) - qJ(4);
t5 = sin(pkin(10));
t6 = cos(pkin(10));
t13 = t5 * t9 - t6 * t7;
t12 = t5 * t7 + t6 * t9;
t11 = qJ(3) * t9 - t16 * t7;
t10 = cos(qJ(1));
t8 = sin(qJ(1));
t4 = t12 * t10;
t3 = t13 * t10;
t2 = t12 * t8;
t1 = t13 * t8;
t14 = [-t2 * r_i_i_C(1) + t1 * r_i_i_C(2) + t15 * t10 - t17 * t8, t3 * r_i_i_C(1) + t4 * r_i_i_C(2) + t11 * t10, t10 * t7, -t8, 0, 0; t4 * r_i_i_C(1) - t3 * r_i_i_C(2) + t17 * t10 + t15 * t8, t1 * r_i_i_C(1) + t2 * r_i_i_C(2) + t11 * t8, t8 * t7, t10, 0, 0; 0 (r_i_i_C(1) * t5 + r_i_i_C(2) * t6 + qJ(3)) * t7 + (r_i_i_C(1) * t6 - r_i_i_C(2) * t5 + t16) * t9, -t9, 0, 0, 0;];
Ja_transl  = t14;
