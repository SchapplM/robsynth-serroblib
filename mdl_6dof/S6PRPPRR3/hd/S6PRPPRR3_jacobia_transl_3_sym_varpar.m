% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% Ja_transl [3x6]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S6PRPPRR3_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobia_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPPRR3_jacobia_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobia_transl_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:59
% EndTime: 2019-02-26 19:45:59
% DurationCPUTime: 0.06s
% Computational Cost: add. (20->10), mult. (47->16), div. (0->0), fcn. (60->6), ass. (0->15)
t17 = pkin(2) + r_i_i_C(1);
t10 = sin(qJ(2));
t6 = sin(pkin(10));
t16 = t6 * t10;
t11 = cos(qJ(2));
t15 = t6 * t11;
t8 = cos(pkin(10));
t14 = t8 * t10;
t13 = t8 * t11;
t12 = r_i_i_C(3) + qJ(3);
t9 = cos(pkin(6));
t7 = sin(pkin(6));
t3 = t9 * t15 + t14;
t1 = -t9 * t13 + t16;
t2 = [0, t12 * (-t9 * t16 + t13) - t17 * t3, t3, 0, 0, 0; 0, t12 * (t9 * t14 + t15) - t17 * t1, t1, 0, 0, 0; 1 (t12 * t10 + t17 * t11) * t7, -t7 * t11, 0, 0, 0;];
Ja_transl  = t2;
