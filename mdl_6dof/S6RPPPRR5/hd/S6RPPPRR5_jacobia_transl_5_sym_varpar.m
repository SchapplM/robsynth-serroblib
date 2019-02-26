% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPPRR5_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR5_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_jacobia_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:25:03
% EndTime: 2019-02-26 20:25:03
% DurationCPUTime: 0.06s
% Computational Cost: add. (34->16), mult. (60->18), div. (0->0), fcn. (78->6), ass. (0->15)
t15 = pkin(7) + r_i_i_C(3);
t14 = pkin(1) + qJ(3);
t13 = pkin(3) + qJ(2);
t12 = sin(pkin(9));
t5 = sin(qJ(5));
t7 = cos(qJ(5));
t11 = t7 * r_i_i_C(1) - t5 * r_i_i_C(2);
t10 = r_i_i_C(1) * t5 + r_i_i_C(2) * t7;
t9 = pkin(4) + t11;
t8 = cos(qJ(1));
t6 = sin(qJ(1));
t4 = cos(pkin(9));
t2 = t8 * t12 + t6 * t4;
t1 = -t6 * t12 + t8 * t4;
t3 = [t9 * t1 + t13 * t8 - t14 * t6 + t15 * t2, t6, t8, 0, -t10 * t2, 0; -t15 * t1 + t13 * t6 + t14 * t8 + t9 * t2, -t8, t6, 0, t10 * t1, 0; 0, 0, 0, 1, t11, 0;];
Ja_transl  = t3;
