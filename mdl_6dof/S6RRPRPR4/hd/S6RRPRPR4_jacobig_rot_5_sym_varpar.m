% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRPR4_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobig_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:39:32
% EndTime: 2019-02-26 21:39:32
% DurationCPUTime: 0.08s
% Computational Cost: add. (8->5), mult. (22->12), div. (0->0), fcn. (35->8), ass. (0->12)
t100 = sin(pkin(11));
t102 = cos(pkin(11));
t104 = sin(qJ(2));
t106 = cos(qJ(2));
t108 = t100 * t104 - t102 * t106;
t107 = cos(qJ(1));
t105 = sin(qJ(1));
t103 = cos(pkin(6));
t101 = sin(pkin(6));
t99 = -t100 * t106 - t102 * t104;
t98 = t108 * t103;
t1 = [0, t105 * t101, 0, -t105 * t98 - t107 * t99, 0, 0; 0, -t107 * t101, 0, -t105 * t99 + t107 * t98, 0, 0; 1, t103, 0, t108 * t101, 0, 0;];
Jg_rot  = t1;