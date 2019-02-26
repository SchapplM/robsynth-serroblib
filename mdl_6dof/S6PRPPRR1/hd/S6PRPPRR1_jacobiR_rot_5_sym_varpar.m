% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPPRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:44:35
% EndTime: 2019-02-26 19:44:35
% DurationCPUTime: 0.08s
% Computational Cost: add. (62->16), mult. (122->36), div. (0->0), fcn. (178->10), ass. (0->24)
t102 = sin(pkin(10));
t103 = sin(pkin(6));
t111 = t102 * t103;
t105 = cos(pkin(10));
t110 = t103 * t105;
t106 = cos(pkin(6));
t101 = sin(pkin(11));
t104 = cos(pkin(11));
t107 = sin(qJ(2));
t108 = cos(qJ(2));
t109 = t108 * t101 + t107 * t104;
t94 = t109 * t106;
t95 = t107 * t101 - t108 * t104;
t88 = -t102 * t95 + t105 * t94;
t90 = -t102 * t94 - t105 * t95;
t100 = pkin(12) + qJ(5);
t99 = cos(t100);
t98 = sin(t100);
t93 = t95 * t106;
t92 = t109 * t103;
t91 = t95 * t103;
t89 = t102 * t93 - t105 * t109;
t87 = -t102 * t109 - t105 * t93;
t1 = [0, t89 * t99, 0, 0, t111 * t99 - t90 * t98, 0; 0, t87 * t99, 0, 0, -t110 * t99 - t88 * t98, 0; 0, -t91 * t99, 0, 0, t106 * t99 - t92 * t98, 0; 0, -t89 * t98, 0, 0, -t111 * t98 - t90 * t99, 0; 0, -t87 * t98, 0, 0, t110 * t98 - t88 * t99, 0; 0, t91 * t98, 0, 0, -t106 * t98 - t92 * t99, 0; 0, t90, 0, 0, 0, 0; 0, t88, 0, 0, 0, 0; 0, t92, 0, 0, 0, 0;];
JR_rot  = t1;
