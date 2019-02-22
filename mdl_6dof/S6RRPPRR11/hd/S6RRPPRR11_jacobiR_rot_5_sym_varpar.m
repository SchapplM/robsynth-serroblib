% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:19
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPPRR11_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:18:58
% EndTime: 2019-02-22 11:18:58
% DurationCPUTime: 0.05s
% Computational Cost: add. (52->16), mult. (87->30), div. (0->0), fcn. (134->8), ass. (0->26)
t92 = sin(pkin(6));
t94 = sin(qJ(2));
t107 = t92 * t94;
t95 = sin(qJ(1));
t106 = t92 * t95;
t96 = cos(qJ(2));
t105 = t92 * t96;
t97 = cos(qJ(1));
t104 = t92 * t97;
t103 = t95 * t94;
t102 = t95 * t96;
t101 = t97 * t94;
t100 = t97 * t96;
t93 = cos(pkin(6));
t83 = -t93 * t100 + t103;
t91 = pkin(11) + qJ(5);
t89 = sin(t91);
t90 = cos(t91);
t99 = t90 * t104 - t83 * t89;
t98 = t89 * t104 + t83 * t90;
t86 = -t93 * t103 + t100;
t85 = t93 * t102 + t101;
t84 = t93 * t101 + t102;
t82 = t90 * t106 + t85 * t89;
t81 = -t89 * t106 + t85 * t90;
t1 = [t99, t86 * t89, 0, 0, t81, 0; t82, t84 * t89, 0, 0, t98, 0; 0, t89 * t107, 0, 0, -t90 * t105 - t93 * t89, 0; -t98, t86 * t90, 0, 0, -t82, 0; t81, t84 * t90, 0, 0, t99, 0; 0, t90 * t107, 0, 0, t89 * t105 - t93 * t90, 0; -t84, -t85, 0, 0, 0, 0; t86, -t83, 0, 0, 0, 0; 0, t105, 0, 0, 0, 0;];
JR_rot  = t1;
