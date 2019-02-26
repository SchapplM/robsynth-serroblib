% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP12_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:52:13
% EndTime: 2019-02-26 21:52:13
% DurationCPUTime: 0.04s
% Computational Cost: add. (56->19), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->18)
t102 = cos(qJ(2));
t99 = qJ(4) + qJ(5);
t97 = sin(t99);
t109 = t102 * t97;
t98 = cos(t99);
t108 = t102 * t98;
t100 = sin(qJ(2));
t101 = sin(qJ(1));
t107 = t101 * t100;
t106 = t101 * t102;
t103 = cos(qJ(1));
t105 = t103 * t100;
t104 = t103 * t102;
t94 = t103 * t98 - t97 * t107;
t93 = t103 * t97 + t98 * t107;
t92 = t101 * t98 + t97 * t105;
t91 = t101 * t97 - t98 * t105;
t1 = [t94, t97 * t104, 0, -t91, -t91, 0; t92, t97 * t106, 0, t93, t93, 0; 0, t100 * t97, 0, -t108, -t108, 0; -t106, -t105, 0, 0, 0, 0; t104, -t107, 0, 0, 0, 0; 0, t102, 0, 0, 0, 0; t93, -t98 * t104, 0, t92, t92, 0; t91, -t98 * t106, 0, -t94, -t94, 0; 0, -t100 * t98, 0, -t109, -t109, 0;];
JR_rot  = t1;
