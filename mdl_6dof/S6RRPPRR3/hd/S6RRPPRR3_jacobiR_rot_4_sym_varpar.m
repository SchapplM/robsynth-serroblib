% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR3_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_jacobiR_rot_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:29:33
% EndTime: 2019-02-26 21:29:33
% DurationCPUTime: 0.05s
% Computational Cost: add. (46->15), mult. (126->32), div. (0->0), fcn. (184->10), ass. (0->22)
t100 = sin(pkin(6));
t105 = sin(qJ(1));
t111 = t100 * t105;
t107 = cos(qJ(1));
t110 = t100 * t107;
t102 = cos(pkin(11));
t104 = sin(qJ(2));
t106 = cos(qJ(2));
t99 = sin(pkin(11));
t109 = t106 * t102 - t104 * t99;
t103 = cos(pkin(6));
t108 = t104 * t102 + t106 * t99;
t94 = t108 * t103;
t89 = -t105 * t109 - t107 * t94;
t91 = -t105 * t94 + t107 * t109;
t101 = cos(pkin(12));
t98 = sin(pkin(12));
t93 = t109 * t103;
t92 = t109 * t100;
t90 = -t105 * t93 - t107 * t108;
t88 = -t105 * t108 + t107 * t93;
t1 = [t89 * t101 + t98 * t110, t90 * t101, 0, 0, 0, 0; t91 * t101 + t98 * t111, t88 * t101, 0, 0, 0, 0; 0, t92 * t101, 0, 0, 0, 0; t101 * t110 - t89 * t98, -t90 * t98, 0, 0, 0, 0; t101 * t111 - t91 * t98, -t88 * t98, 0, 0, 0, 0; 0, -t92 * t98, 0, 0, 0, 0; t88, t91, 0, 0, 0, 0; -t90, -t89, 0, 0, 0, 0; 0, t108 * t100, 0, 0, 0, 0;];
JR_rot  = t1;
