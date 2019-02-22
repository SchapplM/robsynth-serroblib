% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:29
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPRPRP8_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_jacobiR_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:29:16
% EndTime: 2019-02-22 10:29:16
% DurationCPUTime: 0.03s
% Computational Cost: add. (37->15), mult. (40->20), div. (0->0), fcn. (69->6), ass. (0->16)
t87 = sin(qJ(5));
t88 = sin(qJ(1));
t94 = t88 * t87;
t89 = cos(qJ(5));
t93 = t88 * t89;
t90 = cos(qJ(1));
t92 = t90 * t87;
t91 = t90 * t89;
t86 = qJ(3) + pkin(9);
t85 = cos(t86);
t84 = sin(t86);
t83 = t84 * t91 - t94;
t82 = t84 * t92 + t93;
t81 = t84 * t93 + t92;
t80 = t84 * t94 - t91;
t1 = [t83, 0, t85 * t93, 0, -t80, 0; t81, 0, -t85 * t91, 0, t82, 0; 0, 0, -t84 * t89, 0, -t85 * t87, 0; -t90 * t85, 0, t88 * t84, 0, 0, 0; -t88 * t85, 0, -t90 * t84, 0, 0, 0; 0, 0, t85, 0, 0, 0; t82, 0, t85 * t94, 0, t81, 0; t80, 0, -t85 * t92, 0, -t83, 0; 0, 0, -t84 * t87, 0, t85 * t89, 0;];
JR_rot  = t1;
