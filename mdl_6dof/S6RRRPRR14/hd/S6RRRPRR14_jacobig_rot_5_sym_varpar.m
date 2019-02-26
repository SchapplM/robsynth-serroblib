% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR14_jacobig_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_jacobig_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_jacobig_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:45
% EndTime: 2019-02-26 22:23:45
% DurationCPUTime: 0.02s
% Computational Cost: add. (9->9), mult. (24->19), div. (0->0), fcn. (40->8), ass. (0->14)
t100 = sin(qJ(3));
t98 = sin(pkin(6));
t110 = t98 * t100;
t101 = sin(qJ(2));
t102 = sin(qJ(1));
t109 = t102 * t101;
t104 = cos(qJ(2));
t108 = t102 * t104;
t105 = cos(qJ(1));
t107 = t105 * t101;
t106 = t105 * t104;
t103 = cos(qJ(3));
t99 = cos(pkin(6));
t1 = [0, t102 * t98, t99 * t108 + t107, 0 (-t99 * t109 + t106) * t103 + t102 * t110, 0; 0, -t105 * t98, -t99 * t106 + t109, 0 (t99 * t107 + t108) * t103 - t105 * t110, 0; 1, t99, -t98 * t104, 0, t98 * t101 * t103 + t99 * t100, 0;];
Jg_rot  = t1;
