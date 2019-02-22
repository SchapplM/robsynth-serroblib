% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:41
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:40:51
% EndTime: 2019-02-22 11:40:51
% DurationCPUTime: 0.09s
% Computational Cost: add. (147->22), mult. (68->20), div. (0->0), fcn. (117->6), ass. (0->19)
t101 = qJ(4) + qJ(5) + qJ(6);
t97 = sin(t101);
t102 = qJ(2) + pkin(11);
t99 = sin(t102);
t110 = t99 * t97;
t98 = cos(t101);
t109 = t99 * t98;
t103 = sin(qJ(1));
t108 = t103 * t99;
t104 = cos(qJ(1));
t107 = t104 * t99;
t100 = cos(t102);
t106 = t103 * t100;
t105 = t104 * t100;
t96 = t103 * t97 + t105 * t98;
t95 = t103 * t98 - t105 * t97;
t94 = t104 * t97 - t106 * t98;
t93 = t104 * t98 + t106 * t97;
t1 = [t94, -t98 * t107, 0, t95, t95, t95; t96, -t98 * t108, 0, -t93, -t93, -t93; 0, t100 * t98, 0, -t110, -t110, -t110; t93, t97 * t107, 0, -t96, -t96, -t96; t95, t97 * t108, 0, t94, t94, t94; 0, -t100 * t97, 0, -t109, -t109, -t109; -t108, t105, 0, 0, 0, 0; t107, t106, 0, 0, 0, 0; 0, t99, 0, 0, 0, 0;];
JR_rot  = t1;
