% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S4PRRP1
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
%
% Output:
% Ja_transl [3x4]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_transl = S4PRRP1_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_jacobia_transl_3_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRP1_jacobia_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_jacobia_transl_3_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:27:32
% EndTime: 2019-02-26 19:27:32
% DurationCPUTime: 0.05s
% Computational Cost: add. (24->6), mult. (10->6), div. (0->0), fcn. (10->4), ass. (0->7)
t5 = pkin(6) + qJ(2);
t4 = qJ(3) + t5;
t2 = sin(t4);
t3 = cos(t4);
t7 = r_i_i_C(1) * t3 - r_i_i_C(2) * t2;
t6 = -r_i_i_C(1) * t2 - r_i_i_C(2) * t3;
t1 = [0, -pkin(2) * sin(t5) + t6, t6, 0; 0, pkin(2) * cos(t5) + t7, t7, 0; 1, 0, 0, 0;];
Ja_transl  = t1;
