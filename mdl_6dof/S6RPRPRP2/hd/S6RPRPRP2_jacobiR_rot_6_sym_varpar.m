% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:26
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:26:30
% EndTime: 2019-02-22 10:26:30
% DurationCPUTime: 0.04s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t93 = qJ(1) + pkin(9);
t89 = sin(t93);
t94 = sin(qJ(5));
t101 = t89 * t94;
t95 = cos(qJ(5));
t100 = t89 * t95;
t92 = qJ(3) + pkin(10);
t90 = cos(t92);
t99 = t90 * t94;
t98 = t90 * t95;
t91 = cos(t93);
t97 = t91 * t94;
t96 = t91 * t95;
t88 = sin(t92);
t87 = t90 * t96 + t101;
t86 = t90 * t97 - t100;
t85 = t89 * t98 - t97;
t84 = -t89 * t99 - t96;
t1 = [-t85, 0, -t88 * t96, 0, -t86, 0; t87, 0, -t88 * t100, 0, t84, 0; 0, 0, t98, 0, -t88 * t94, 0; -t89 * t88, 0, t91 * t90, 0, 0, 0; t91 * t88, 0, t89 * t90, 0, 0, 0; 0, 0, t88, 0, 0, 0; t84, 0, -t88 * t97, 0, t87, 0; t86, 0, -t88 * t101, 0, t85, 0; 0, 0, t99, 0, t88 * t95, 0;];
JR_rot  = t1;
