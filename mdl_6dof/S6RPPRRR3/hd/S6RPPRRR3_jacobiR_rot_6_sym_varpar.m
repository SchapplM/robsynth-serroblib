% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:18
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRR3_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_jacobiR_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:18:47
% EndTime: 2019-02-22 10:18:47
% DurationCPUTime: 0.04s
% Computational Cost: add. (88->19), mult. (54->20), div. (0->0), fcn. (93->6), ass. (0->17)
t87 = qJ(5) + qJ(6);
t84 = sin(t87);
t88 = sin(qJ(4));
t93 = t88 * t84;
t85 = cos(t87);
t92 = t88 * t85;
t89 = cos(qJ(4));
t91 = t89 * t84;
t90 = t89 * t85;
t86 = qJ(1) + pkin(10);
t83 = cos(t86);
t82 = sin(t86);
t81 = -t82 * t84 + t83 * t92;
t80 = t82 * t85 + t83 * t93;
t79 = t82 * t92 + t83 * t84;
t78 = -t82 * t93 + t83 * t85;
t1 = [t81, 0, 0, t82 * t90, t78, t78; t79, 0, 0, -t83 * t90, t80, t80; 0, 0, 0, -t92, -t91, -t91; -t80, 0, 0, -t82 * t91, -t79, -t79; t78, 0, 0, t83 * t91, t81, t81; 0, 0, 0, t93, -t90, -t90; -t83 * t89, 0, 0, t82 * t88, 0, 0; -t82 * t89, 0, 0, -t83 * t88, 0, 0; 0, 0, 0, t89, 0, 0;];
JR_rot  = t1;
