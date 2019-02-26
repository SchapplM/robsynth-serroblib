% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RPPR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
%
% Output:
% Ja_transl [3x4]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4RPPR1_jacobia_transl_4_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_jacobia_transl_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPR1_jacobia_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_jacobia_transl_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:30:33
% EndTime: 2019-02-26 19:30:33
% DurationCPUTime: 0.06s
% Computational Cost: add. (45->13), mult. (32->14), div. (0->0), fcn. (42->6), ass. (0->11)
t13 = pkin(2) + pkin(3);
t12 = cos(qJ(4));
t11 = sin(qJ(4));
t8 = qJ(1) + pkin(6);
t6 = sin(t8);
t7 = cos(t8);
t1 = -t6 * t11 - t7 * t12;
t2 = t7 * t11 - t6 * t12;
t10 = -t2 * r_i_i_C(1) + t1 * r_i_i_C(2);
t9 = t1 * r_i_i_C(1) + t2 * r_i_i_C(2);
t3 = [-sin(qJ(1)) * pkin(1) + t7 * qJ(3) - t13 * t6 - t10, 0, t6, t10; cos(qJ(1)) * pkin(1) + t6 * qJ(3) + t13 * t7 - t9, 0, -t7, t9; 0, 1, 0, 0;];
Ja_transl  = t3;
