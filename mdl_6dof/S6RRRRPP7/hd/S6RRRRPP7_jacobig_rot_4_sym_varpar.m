% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPP7_jacobig_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_jacobig_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_jacobig_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:29:04
% EndTime: 2019-02-26 22:29:04
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->15)
t101 = sin(pkin(6));
t105 = sin(qJ(1));
t114 = t105 * t101;
t104 = sin(qJ(2));
t113 = t105 * t104;
t107 = cos(qJ(2));
t112 = t105 * t107;
t108 = cos(qJ(1));
t111 = t108 * t101;
t110 = t108 * t104;
t109 = t108 * t107;
t106 = cos(qJ(3));
t103 = sin(qJ(3));
t102 = cos(pkin(6));
t1 = [0, t114, t102 * t112 + t110 (-t102 * t113 + t109) * t103 - t106 * t114, 0, 0; 0, -t111, -t102 * t109 + t113 (t102 * t110 + t112) * t103 + t106 * t111, 0, 0; 1, t102, -t101 * t107, t101 * t104 * t103 - t102 * t106, 0, 0;];
Jg_rot  = t1;
