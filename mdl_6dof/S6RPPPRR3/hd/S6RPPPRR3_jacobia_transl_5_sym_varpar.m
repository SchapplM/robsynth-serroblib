% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6RPPPRR3_jacobia_transl_5_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobia_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR3_jacobia_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobia_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:23:51
% EndTime: 2019-02-26 20:23:51
% DurationCPUTime: 0.07s
% Computational Cost: add. (48->17), mult. (64->19), div. (0->0), fcn. (84->7), ass. (0->15)
t19 = pkin(1) + pkin(2);
t18 = r_i_i_C(3) + pkin(7) + qJ(4);
t17 = cos(qJ(1));
t16 = sin(qJ(1));
t15 = cos(pkin(9));
t14 = sin(pkin(9));
t9 = pkin(10) + qJ(5);
t7 = sin(t9);
t8 = cos(t9);
t13 = -t8 * r_i_i_C(1) + t7 * r_i_i_C(2);
t12 = r_i_i_C(1) * t7 + r_i_i_C(2) * t8;
t11 = -t13 + cos(pkin(10)) * pkin(4) + pkin(3);
t2 = t17 * t14 - t16 * t15;
t1 = -t16 * t14 - t17 * t15;
t3 = [t17 * qJ(2) + t18 * t1 + t11 * t2 - t19 * t16, t16, 0, t2, t12 * t1, 0; t16 * qJ(2) - t11 * t1 + t19 * t17 + t18 * t2, -t17, 0, -t1, t12 * t2, 0; 0, 0, -1, 0, t13, 0;];
Ja_transl  = t3;
