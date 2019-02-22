% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:14
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:13:55
% EndTime: 2019-02-22 10:13:55
% DurationCPUTime: 0.03s
% Computational Cost: add. (59->14), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->19)
t90 = qJ(1) + pkin(9);
t86 = sin(t90);
t91 = sin(qJ(5));
t98 = t86 * t91;
t92 = cos(qJ(5));
t97 = t86 * t92;
t89 = pkin(10) + qJ(4);
t87 = cos(t89);
t96 = t87 * t91;
t95 = t87 * t92;
t88 = cos(t90);
t94 = t88 * t91;
t93 = t88 * t92;
t85 = sin(t89);
t84 = t87 * t93 + t98;
t83 = t87 * t94 - t97;
t82 = t86 * t95 - t94;
t81 = -t86 * t96 - t93;
t1 = [-t82, 0, 0, -t85 * t93, -t83, 0; t84, 0, 0, -t85 * t97, t81, 0; 0, 0, 0, t95, -t85 * t91, 0; -t86 * t85, 0, 0, t88 * t87, 0, 0; t88 * t85, 0, 0, t86 * t87, 0, 0; 0, 0, 0, t85, 0, 0; t81, 0, 0, -t85 * t94, t84, 0; t83, 0, 0, -t85 * t98, t82, 0; 0, 0, 0, t96, t85 * t92, 0;];
JR_rot  = t1;
