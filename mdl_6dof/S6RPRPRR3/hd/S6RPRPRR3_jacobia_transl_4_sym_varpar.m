% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR3
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
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPRPRR3_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_jacobia_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR3_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_jacobia_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:19
% EndTime: 2019-02-26 20:50:19
% DurationCPUTime: 0.07s
% Computational Cost: add. (61->14), mult. (67->18), div. (0->0), fcn. (74->8), ass. (0->14)
t4 = sin(pkin(11));
t5 = cos(pkin(11));
t11 = r_i_i_C(1) * t5 - r_i_i_C(2) * t4 + pkin(3);
t12 = r_i_i_C(3) + qJ(4);
t6 = sin(qJ(3));
t7 = cos(qJ(3));
t9 = t11 * t7 + t12 * t6;
t13 = pkin(2) + t9;
t10 = r_i_i_C(1) * t4 + r_i_i_C(2) * t5 + pkin(7);
t8 = -t11 * t6 + t12 * t7;
t3 = qJ(1) + pkin(10);
t2 = cos(t3);
t1 = sin(t3);
t14 = [-sin(qJ(1)) * pkin(1) + t10 * t2 - t13 * t1, 0, t8 * t2, t2 * t6, 0, 0; cos(qJ(1)) * pkin(1) + t10 * t1 + t13 * t2, 0, t8 * t1, t1 * t6, 0, 0; 0, 1, t9, -t7, 0, 0;];
Ja_transl  = t14;
