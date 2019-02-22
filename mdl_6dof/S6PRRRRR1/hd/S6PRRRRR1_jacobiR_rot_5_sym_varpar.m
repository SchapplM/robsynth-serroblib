% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:02
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:02:00
% EndTime: 2019-02-22 10:02:00
% DurationCPUTime: 0.05s
% Computational Cost: add. (123->14), mult. (117->32), div. (0->0), fcn. (180->8), ass. (0->26)
t97 = sin(pkin(12));
t98 = sin(pkin(6));
t108 = t97 * t98;
t99 = cos(pkin(12));
t107 = t98 * t99;
t102 = cos(qJ(2));
t106 = t102 * t98;
t101 = sin(qJ(2));
t105 = t98 * t101;
t100 = cos(pkin(6));
t104 = t100 * t101;
t103 = t100 * t102;
t96 = qJ(3) + qJ(4) + qJ(5);
t95 = cos(t96);
t94 = sin(t96);
t93 = t99 * t102 - t97 * t104;
t92 = -t99 * t101 - t97 * t103;
t91 = t97 * t102 + t99 * t104;
t90 = -t97 * t101 + t99 * t103;
t89 = -t100 * t94 - t95 * t105;
t88 = t100 * t95 - t94 * t105;
t87 = -t94 * t108 - t93 * t95;
t86 = t95 * t108 - t93 * t94;
t85 = t94 * t107 - t91 * t95;
t84 = -t95 * t107 - t91 * t94;
t1 = [0, t92 * t95, t86, t86, t86, 0; 0, t90 * t95, t84, t84, t84, 0; 0, t95 * t106, t88, t88, t88, 0; 0, -t92 * t94, t87, t87, t87, 0; 0, -t90 * t94, t85, t85, t85, 0; 0, -t94 * t106, t89, t89, t89, 0; 0, t93, 0, 0, 0, 0; 0, t91, 0, 0, 0, 0; 0, t105, 0, 0, 0, 0;];
JR_rot  = t1;
