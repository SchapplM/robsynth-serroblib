% Translatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4PRRP3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
%
% Output:
% Ja_transl [3x4]
%   Translatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:12
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_transl = S4PRRP3_jacobia_transl_4_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_jacobia_transl_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4PRRP3_jacobia_transl_4_floatb_twist_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP3_jacobia_transl_4_floatb_twist_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobia_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:12:40
% EndTime: 2018-11-14 14:12:40
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->6), mult. (14->6), div. (0->0), fcn. (14->4), ass. (0->7)
t8 = pkin(3) + r_i_i_C(1);
t5 = qJ(2) + qJ(3);
t3 = sin(t5);
t4 = cos(t5);
t7 = -t3 * r_i_i_C(2) + t8 * t4;
t6 = -t4 * r_i_i_C(2) - t8 * t3;
t1 = [1, cos(qJ(2)) * pkin(2) + t7, t7, 0; 0, -sin(qJ(2)) * pkin(2) + t6, t6, 0; 0, 0, 0, -1;];
Ja_transl  = t1;
